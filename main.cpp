#include <iostream>
#include "pointmatcher/PointMatcher.h"

typedef PointMatcher<float> PM;

static double gicp_epsilon = 1e-3;

int main(int argc, char** argv)
{
	PM::DataPoints reading = PM::DataPoints::load("/home/norlab/Desktop/1617238752456772327.vtk");
	PM::DataPoints reference = PM::DataPoints::load("/home/norlab/Desktop/1617238750956836700.vtk");

	PM::ICP icp;
	icp.setDefault();
	std::shared_ptr<PM::ErrorMinimizer> errorMinimizer = PM::get().ErrorMinimizerRegistrar.create("GICPErrorMinimizer");
	icp.errorMinimizer = errorMinimizer;
	std::cout << icp(reading, reference) << std::endl;

	return 0;
}
