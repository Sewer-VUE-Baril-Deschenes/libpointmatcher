#include <iostream>
#include "pointmatcher/PointMatcher.h"

typedef PointMatcher<float> PM;

int main(int argc, char** argv)
{
	PM::DataPoints reading = PM::DataPoints::load("/home/norlab/Desktop/1617238752456772327.vtk");
	PM::DataPoints reference = PM::DataPoints::load("/home/norlab/Desktop/1617238750956836700.vtk");

	std::shared_ptr<PM::Transformation> transformator = PM::get().TransformationRegistrar.create("RigidTransformation");
	PM::TransformationParameters perturbation = PM::TransformationParameters::Identity(4, 4);
	perturbation(2, 3) = 0.3;
	PM::DataPoints perturbatedReading = transformator->compute(reading, perturbation);
	perturbatedReading.save("/home/norlab/Desktop/perturbated_1617238752456772327.vtk");

	std::shared_ptr<PM::DataPointsFilter> gicpCovarianceFilter = PM::get().DataPointsFilterRegistrar.create("SurfaceNormalGICPCovarianceDataPointsFilter");
	gicpCovarianceFilter->inPlaceFilter(reference);

	PM::ICP icp;
	icp.setDefault();
	std::shared_ptr<PM::ErrorMinimizer> errorMinimizer = PM::get().ErrorMinimizerRegistrar.create("GICPErrorMinimizer");
	icp.errorMinimizer = errorMinimizer;

	PM::TransformationParameters correction = icp(perturbatedReading, reference);
	std::cout << correction << std::endl;

	PM::DataPoints transformedReading = transformator->compute(perturbatedReading, correction);
	transformedReading.save("/home/norlab/Desktop/registered_1617238752456772327.vtk");

	return 0;
}
