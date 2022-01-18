#include <iostream>
#include "pointmatcher/PointMatcher.h"
#include "pointmatcher/ErrorMinimizers/GICP/GICP.h"

typedef PointMatcher<float> PM;

static double gicp_epsilon = 1e-3;

int main(int argc, char** argv)
{
	PM::DataPoints reading = PM::DataPoints::load("/home/norlab/Desktop/1617238752456772327.vtk");
	PM::DataPoints reference = PM::DataPoints::load("/home/norlab/Desktop/1617238750956836700.vtk");

	std::cout << reference.descriptorLabels << std::endl;

	return 0;
}
