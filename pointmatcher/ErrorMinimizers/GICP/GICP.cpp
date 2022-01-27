/*************************************************************
  Generalized-ICP Copyright (c) 2009 Aleksandr Segal.
  All rights reserved.

  Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that the
  following conditions are met:

* Redistributions of source code must retain the above
  copyright notice, this list of conditions and the
  following disclaimer.
* Redistributions in binary form must reproduce the above
  copyright notice, this list of conditions and the
  following disclaimer in the documentation and/or other
  materials provided with the distribution.
* The names of the contributors may not be used to endorse
  or promote products derived from this software
  without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.
*************************************************************/

#include "GICP.h"
#include "transform.h"
#include "optimize.h"

template<typename T>
GICPErrorMinimizer<T>::GICPErrorMinimizer(const Parameters& params)
{
	// initialize some stuff
}

template<typename T>
GICPErrorMinimizer<T>::GICPErrorMinimizer(const ParametersDoc paramsDoc, const Parameters& params)
{
	// initialize some stuff
}


template<typename T>
void convertDescriptorToGSLMatrix(const typename PointMatcher<T>::Vector& descriptor, gsl_matrix* GSLMatrix)
{
	for(unsigned int i = 0; i < GSLMatrix->size1; ++i)
	{
		for(unsigned int j = 0; j < GSLMatrix->size2; ++j)
		{
			gsl_matrix_set(GSLMatrix, i, j, descriptor(j * GSLMatrix->size1 + i));
		}
	}
}

template<typename T>
typename PointMatcher<T>::TransformationParameters convertDGCTransformToTransformationParameters(const dgc_transform_t& DGCTransform, unsigned int nbRows, unsigned int nbCols)
{
	typename PointMatcher<T>::TransformationParameters transformationParameters(nbRows, nbCols);
	for(unsigned int i = 0; i < nbRows; ++i)
	{
		for(unsigned int j = 0; j < nbCols; ++j)
		{
			transformationParameters(i, j) = DGCTransform[i][j];
		}
	}
	return transformationParameters;
}

template<typename T>
typename PointMatcher<T>::TransformationParameters GICPErrorMinimizer<T>::compute(const ErrorElements& mPtsIn)
{
	if (!mPtsIn.reading.descriptorExists("gicpCovariance"))
	{
		throw InvalidField("GICPErrorMinimizer: Error, no gicpCovariance found in reading descriptors.");
	}
	if (!mPtsIn.reference.descriptorExists("gicpCovariance"))
	{
		throw InvalidField("GICPErrorMinimizer: Error, no gicpCovariance found in reference descriptors.");
	}
	ErrorElements mPts = mPtsIn;
	const auto& readingCovariances = mPts.reading.getDescriptorViewByName("gicpCovariance");
	const auto& referenceCovariances = mPts.reference.getDescriptorViewByName("gicpCovariance");

	int n = mPts.reading.getNbPoints();
	dgc_transform_t t;
	dgc_transform_identity(t);
	dgc_transform_t identity;
	dgc_transform_identity(identity);
	double query_point[3];

	dgc::gicp::gicp_mat_t* mahalanobis = new dgc::gicp::gicp_mat_t[n];

	gsl_matrix* gsl_R = gsl_matrix_alloc(3, 3);
	gsl_matrix* gsl_temp = gsl_matrix_alloc(3, 3);
	gsl_matrix* C1 = gsl_matrix_alloc(3, 3);
	gsl_matrix* C2 = gsl_matrix_alloc(3, 3);

	/* set up the optimization parameters */
	dgc::gicp::GICPOptData<T> opt_data;
	opt_data.nn_indecies = &mPts.matches.ids;
	opt_data.p1 = &mPts.reading;
	opt_data.p2 = &mPts.reference;
	opt_data.M = mahalanobis;
	opt_data.solve_rotation = true;
	opt_data.num_matches = n;
	dgc_transform_copy(opt_data.base_t, identity);

	dgc::gicp::GICPOptimizer<T> opt;
	opt.SetDebug(false);
	opt.SetMaxIterations(10);
	/* set up the mahalanobis matricies */
	/* these are identity for now to ease debugging */
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < 3; k++)
		{
			for(int l = 0; l < 3; l++)
			{
				mahalanobis[i][k][l] = (k == l) ? 1 : 0.;
			}
		}
	}

	dgc_transform_t transform_R;
	dgc_transform_copy(transform_R, identity);
	dgc_transform_left_multiply(transform_R, t);
	// copy the rotation component of the current total transformation (including base), into a gsl matrix
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			gsl_matrix_set(gsl_R, i, j, transform_R[i][j]);
		}
	}

	for(int i = 0; i < n; i++)
	{
		query_point[0] = mPts.reading.features(0, i);
		query_point[1] = mPts.reading.features(1, i);
		query_point[2] = mPts.reading.features(2, i);

		dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], identity);
		dgc_transform_point(&query_point[0], &query_point[1], &query_point[2], t);

		// set up the updated mahalanobis matrix here
		convertDescriptorToGSLMatrix<T>(readingCovariances.col(i), C1); // Plane-to-plane
//		gsl_matrix_set_zero(C1); // Point-to-plane
		convertDescriptorToGSLMatrix<T>(referenceCovariances.col(i), C2);
		gsl_matrix_view M = gsl_matrix_view_array(&mahalanobis[i][0][0], 3, 3);
		gsl_matrix_set_zero(&M.matrix);
		gsl_matrix_set_zero(gsl_temp);

		// M = R*C1  // using M as a temp variable here
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., gsl_R, C1, 1., &M.matrix);

		// temp = M*R' // move the temp value to 'temp' here
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., &M.matrix, gsl_R, 0., gsl_temp);

		// temp += C2
		gsl_matrix_add(gsl_temp, C2);
		// at this point temp = C2 + R*C1*R'

		// now invert temp to get the mahalanobis distance metric for gicp
		// M = temp^-1
		gsl_matrix_set_identity(&M.matrix);
		gsl_error_handler_t* error_handler = gsl_set_error_handler_off();
		int status = gsl_linalg_cholesky_decomp(gsl_temp);
		if(status != GSL_EDOM)
		{
			for(int k = 0; k < 3; k++)
			{
				gsl_vector_view row_view = gsl_matrix_row(&M.matrix, k);
				gsl_linalg_cholesky_svx(gsl_temp, &row_view.vector);
			}
		}
		gsl_set_error_handler(error_handler);
	}

	/* optimize transformation using the current assignment and Mahalanobis metrics*/
	opt.Optimize(t, opt_data);

	if(mahalanobis != NULL)
	{
		delete[] mahalanobis;
	}
	if(gsl_R != NULL)
	{
		gsl_matrix_free(gsl_R);
	}
	if(gsl_temp != NULL)
	{
		gsl_matrix_free(gsl_temp);
	}
	if(C1 != NULL)
	{
		gsl_matrix_free(C1);
	}
	if(C2 != NULL)
	{
		gsl_matrix_free(C2);
	}

	return convertDGCTransformToTransformationParameters<T>(t, 4, 4);
}

template struct GICPErrorMinimizer<float>;
template struct GICPErrorMinimizer<double>;