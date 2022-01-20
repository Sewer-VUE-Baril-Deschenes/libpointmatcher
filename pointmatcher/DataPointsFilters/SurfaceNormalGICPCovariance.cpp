// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2018,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include "SurfaceNormalGICPCovariance.h"

// Eigenvalues
#include "Eigen/QR"
#include "Eigen/Eigenvalues"

#include "PointMatcherPrivate.h"
#include "IO.h"
#include "MatchersImpl.h"

#include <boost/format.hpp>

#include "DataPointsFilters/utils/utils.h"

// SurfaceNormalGICPCovarianceDataPointsFilter
// Constructor
template<typename T>
SurfaceNormalGICPCovarianceDataPointsFilter<T>::SurfaceNormalGICPCovarianceDataPointsFilter(const Parameters& params):
        PointMatcher<T>::DataPointsFilter("SurfaceNormalGICPCovarianceDataPointsFilter",
                                          SurfaceNormalGICPCovarianceDataPointsFilter::availableParameters(), params),
        knn(Parametrizable::get<int>("knn")),
        maxDist(Parametrizable::get<T>("maxDist")),
        epsilon(Parametrizable::get<T>("epsilon")),
        keepNormals(Parametrizable::get<bool>("keepNormals")),
        keepDensities(Parametrizable::get<bool>("keepDensities")),
        keepEigenValues(Parametrizable::get<bool>("keepEigenValues")),
        keepEigenVectors(Parametrizable::get<bool>("keepEigenVectors")),
        keepMatchedIds(Parametrizable::get<bool>("keepMatchedIds")),
        keepMeanDist(Parametrizable::get<bool>("keepMeanDist")),
        keepGICPCovariance(Parametrizable::get<bool>("keepGICPCovariance")),
        sortEigen(Parametrizable::get<bool>("sortEigen")),
        smoothNormals(Parametrizable::get<bool>("smoothNormals")),
        measurementCovariance(Parametrizable::get<T>("measurementCovariance")),
        normalVariance(Parametrizable::get<T>("normalVariance"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints
SurfaceNormalGICPCovarianceDataPointsFilter<T>::filter(
        const DataPoints& input)
{
    DataPoints output(input);
    inPlaceFilter(output);
    return output;
}

template<typename T>
typename PointMatcher<T>::Matrix inverseSquareRootDiagonalMatrix(typename PointMatcher<T>::Matrix matrix)
{
    typename PointMatcher<T>::Matrix inverseSquareRootedMatrix = PointMatcher<T>::Matrix::Zero(matrix.rows(), matrix.cols());
    for(size_t i=0; i < inverseSquareRootedMatrix.rows(); ++i)
    {
        inverseSquareRootedMatrix(i, i) = 1.0 / std::sqrt(matrix(i, i));
    }
    return inverseSquareRootedMatrix;
}

// In-place filter
template<typename T>
void SurfaceNormalGICPCovarianceDataPointsFilter<T>::inPlaceFilter(
        DataPoints& cloud)
{
    typedef typename DataPoints::View View;
    typedef typename DataPoints::Label Label;
    typedef typename DataPoints::Labels Labels;
    typedef typename MatchersImpl<T>::KDTreeMatcher KDTreeMatcher;
    typedef typename PointMatcher<T>::Matches Matches;

    using namespace PointMatcherSupport;

    const int pointsCount(cloud.features.cols());
    const int featDim(cloud.features.rows());
    const int descDim(cloud.descriptors.rows());
    const unsigned int labelDim(cloud.descriptorLabels.size());

    // Validate descriptors and labels
    int insertDim(0);
    for(unsigned int i = 0; i < labelDim ; ++i)
        insertDim += cloud.descriptorLabels[i].span;
    if (insertDim != descDim)
        throw InvalidField("SurfaceNormalCovarianceDataPointsFilter: Error, descriptor labels do not match descriptor data");

    // Reserve memory for new descriptors
    const int dimNormals(featDim-1);
    const int dimDensities(1);
    const int dimEigValues(featDim-1);
    const int dimEigVectors((featDim-1)*(featDim-1));
    //const int dimMatchedIds(knn);
    const int dimMeanDist(1);
    const int dimGICPCovariance((featDim-1)*(featDim-1));

    boost::optional<View> normals;
    boost::optional<View> densities;
    boost::optional<View> eigenValues;
    boost::optional<View> eigenVectors;
    boost::optional<View> matchedValues;
    boost::optional<View> matchIds;
    boost::optional<View> meanDists;
    boost::optional<View> gicpCovariance;

    Labels cloudLabels;
    if (keepNormals)
        cloudLabels.push_back(Label("normals", dimNormals));
    if (keepDensities)
        cloudLabels.push_back(Label("densities", dimDensities));
    if (keepEigenValues)
        cloudLabels.push_back(Label("eigValues", dimEigValues));
    if (keepEigenVectors)
        cloudLabels.push_back(Label("eigVectors", dimEigVectors));
    if (keepMatchedIds)
        cloudLabels.push_back(Label("matchedIds", knn));
    if (keepMeanDist)
        cloudLabels.push_back(Label("meanDists", dimMeanDist));
    if (keepGICPCovariance)
        cloudLabels.push_back(Label("gicpCovariance", dimGICPCovariance));

    // Reserve memory
    cloud.allocateDescriptors(cloudLabels);

    if (keepNormals)
        normals = cloud.getDescriptorViewByName("normals");
    if (keepDensities)
        densities = cloud.getDescriptorViewByName("densities");
    if (keepEigenValues)
        eigenValues = cloud.getDescriptorViewByName("eigValues");
    if (keepEigenVectors)
        eigenVectors = cloud.getDescriptorViewByName("eigVectors");
    if (keepMatchedIds)
        matchIds = cloud.getDescriptorViewByName("matchedIds");
    if (keepMeanDist)
        meanDists = cloud.getDescriptorViewByName("meanDists");
    if (keepGICPCovariance)
        gicpCovariance = cloud.getDescriptorViewByName("gicpCovariance");


    using namespace PointMatcherSupport;
    // Build kd-tree
    Parametrizable::Parameters param;
    boost::assign::insert(param) ( "knn", toParam(knn) );
    boost::assign::insert(param) ( "epsilon", toParam(epsilon) );
    boost::assign::insert(param) ( "maxDist", toParam(maxDist) );

    KDTreeMatcher matcher(param);
    matcher.init(cloud);

    Matches matches(typename Matches::Dists(knn, pointsCount), typename Matches::Ids(knn, pointsCount));
    matches = matcher.findClosests(cloud);

    const auto& intensities = cloud.getDescriptorViewByName("intensity");

    // Search for surrounding points and compute descriptors
    int degenerateCount(0);
    for (int i = 0; i < pointsCount; ++i)
    {
        bool isDegenerate = false;
        // Mean of nearest neighbors (NN)
        Matrix d(featDim-1, knn);
        int realKnn = 0;
        for(int j = 0; j < int(knn); ++j)
        {
            if (matches.dists(j,i) != Matches::InvalidDist)
            {
                const int refIndex(matches.ids(j,i));
                d.col(realKnn) = cloud.features.block(0, refIndex, featDim-1, 1);
                ++realKnn;
            }
        }
        d.conservativeResize(Eigen::NoChange, realKnn);

        const Vector mean = d.rowwise().sum() / T(realKnn);
        const Matrix NN = d.colwise() - mean;

        const Matrix C((NN * NN.transpose()) / T(realKnn));
        Vector eigenVa = Vector::Zero(featDim-1, 1);
        Matrix eigenVe = Matrix::Zero(featDim-1, featDim-1);
        // Ensure that the matrix is suited for eigenvalues calculation
        Matrix pointWeightedCovarianceInScanFrame;
        if(keepGICPCovariance)
        {
            if(C.fullPivHouseholderQr().rank()+1 >= featDim-1) {
                const Eigen::EigenSolver<Matrix> solver(C);
                eigenVa = solver.eigenvalues().real();
                eigenVe = solver.eigenvectors().real();
                pointWeightedCovarianceInScanFrame = Matrix::Zero(featDim - 1, featDim - 1);
                // sort eigen vectors first
                const std::vector<size_t> idx = sortIndexes<T>(eigenVa);
                const size_t idxSize = idx.size();
                Vector tmp_eigenVa = eigenVa;
                Matrix tmp_eigenVe = eigenVe;
                for (size_t j = 0; j < idxSize; ++j) {
                    eigenVa(j, 0) = tmp_eigenVa(idx[j], 0);
                    eigenVe.col(j) = tmp_eigenVe.col(idx[j]);
                }
                Matrix projectionMatrix(2, 3);
                projectionMatrix.row(0) = eigenVe.col(1).transpose();
                projectionMatrix.row(1) = eigenVe.col(2).transpose();
                Matrix descriptorWeights(1, realKnn);
                for (int j = 0; j < int(knn); ++j) {
                    if (matches.ids(0, i) != Matches::InvalidId) {
                        descriptorWeights(0, j) = exp(
                                -0.5 * std::pow((intensities(0, matches.ids(j)) - intensities(0, i)), 2) /
                                measurementCovariance);
                    }
                }
                Matrix projectedNeighbors(2, realKnn);
                Matrix weightedProjectedNeighbors(2, realKnn);
                for (int j = 0; j < int(knn); ++j) {
                    if (matches.ids(0, i) != Matches::InvalidId) {
                        projectedNeighbors.col(j) = projectionMatrix * cloud.features.col(matches.ids(j)).head(3);
                        weightedProjectedNeighbors.col(j) = descriptorWeights(0, j) * projectedNeighbors.col(j);
                    }
                }
                T weightedSumInverse = 1 / descriptorWeights.sum();
                Vector weightedMean = weightedSumInverse * weightedProjectedNeighbors.rowwise().sum();
                Matrix weightedCovariance = Matrix::Zero(featDim - 2, featDim - 2);
                for (int j = 0; j < int(knn); ++j) {
                    if (matches.ids(0, i) != Matches::InvalidId) {
                        Vector projectedNeighborError = projectedNeighbors.col(j) - weightedMean;
                        weightedCovariance +=
                                descriptorWeights(0, j) * (projectedNeighborError * projectedNeighborError.transpose());
                    }
                }
                weightedCovariance = weightedSumInverse * weightedCovariance;
                Matrix populationCovariance = Matrix::Zero(featDim - 2, featDim - 2);
                populationCovariance(0, 0) = eigenVa(1);
                populationCovariance(1, 1) = eigenVa(2);
                Matrix populationCovarianceInverseSquareRoot = inverseSquareRootDiagonalMatrix<T>(populationCovariance);
                weightedCovariance = populationCovarianceInverseSquareRoot * weightedCovariance *
                                     populationCovarianceInverseSquareRoot;
                Matrix pointWeightedCovarianceInPlane = Matrix::Zero(featDim - 1, featDim - 1);
                pointWeightedCovarianceInPlane.topLeftCorner(featDim - 2, featDim - 2) = weightedCovariance;
                pointWeightedCovarianceInPlane(featDim - 2, featDim - 2) = normalVariance;
                Matrix rotationPlaneToScan(featDim - 1, featDim - 1);
                rotationPlaneToScan.col(0) = eigenVe.col(1);
                rotationPlaneToScan.col(1) = eigenVe.col(2);
                rotationPlaneToScan.col(2) = eigenVe.col(
                        0); // NOTE: The eigen vectors order is determined because surface normal is along the Z axis
                pointWeightedCovarianceInScanFrame =
                        rotationPlaneToScan * pointWeightedCovarianceInPlane * rotationPlaneToScan.transpose();
            }
            else
            {
                pointWeightedCovarianceInScanFrame = Matrix::Zero(featDim - 1, featDim - 1);
            }
        }
        if(keepNormals || keepEigenValues || keepEigenVectors)
        {
            if(C.fullPivHouseholderQr().rank()+1 >= featDim-1)
            {
                const Eigen::EigenSolver<Matrix> solver(C);
                eigenVa = solver.eigenvalues().real();
                eigenVe = solver.eigenvectors().real();

                if(sortEigen)
                {
                    const std::vector<size_t> idx = sortIndexes<T>(eigenVa);
                    const size_t idxSize = idx.size();
                    Vector tmp_eigenVa = eigenVa;
                    Matrix tmp_eigenVe = eigenVe;
                    for(size_t i=0; i<idxSize; ++i)
                    {
                        eigenVa(i,0) = tmp_eigenVa(idx[i], 0);
                        eigenVe.col(i) = tmp_eigenVe.col(idx[i]);
                    }
                }
            }
            else
            {
                //std::cout << "WARNING: Matrix C needed for eigen decomposition is degenerated. Expected cause: no noise in data" << std::endl;
                ++degenerateCount;
                isDegenerate = true;
            }
        }

        if(keepNormals)
        {
            if(sortEigen)
                normals->col(i) = eigenVe.col(0);
            else
                normals->col(i) = computeNormal<T>(eigenVa, eigenVe);

            // clamp normals to [-1,1] to handle approximation errors
            normals->col(i) = normals->col(i).cwiseMax(-1.0).cwiseMin(1.0);
        }
        if(keepDensities)
        {
            if(isDegenerate)
                (*densities)(0, i) = 0.;
            else
                (*densities)(0, i) = computeDensity<T>(NN);
        }
        if(keepEigenValues)
            eigenValues->col(i) = eigenVa;
        if(keepEigenVectors)
            eigenVectors->col(i) = serializeEigVec<T>(eigenVe);
        if(keepMeanDist)
        {
            if(isDegenerate)
                (*meanDists)(0, i) = std::numeric_limits<std::size_t>::max();
            else
            {
                const Vector point = cloud.features.block(0, i, featDim-1, 1);
                (*meanDists)(0, i) = (point - mean).norm();
            }
        }
        if(keepGICPCovariance)
        {
            gicpCovariance->col(i) = serializeEigVec<T>(pointWeightedCovarianceInScanFrame);
        }


    }

    if(keepMatchedIds)
    {
        matchIds.get() = matches.ids.template cast<T>();
    }

    if(smoothNormals)
    {
        for (int i = 0; i < pointsCount; ++i)
        {
            const Vector currentNormal = normals->col(i);
            Vector mean = Vector::Zero(featDim-1);
            int n=0;
            for(int j = 0; j < int(knn); ++j)
            {
                if (matches.dists(j,i) != Matches::InvalidDist)
                {
                    const int refIndex(matches.ids(j,i));
                    const Vector normal = normals->col(refIndex);
                    if(currentNormal.dot(normal) > 0.)
                        mean += normal;
                    else // flip normal vector
                        mean -= normal;

                    ++n;
                }
            }

            normals->col(i) = mean / T(n);
        }
    }

    if (degenerateCount)
    {
        LOG_WARNING_STREAM("WARNING: Matrix C needed for eigen decomposition was degenerated in " << degenerateCount << " points over " << pointsCount << " (" << float(degenerateCount)*100.f/float(pointsCount) << " %)");
    }

}

template struct SurfaceNormalGICPCovarianceDataPointsFilter<float>;
template struct SurfaceNormalGICPCovarianceDataPointsFilter<double>;

