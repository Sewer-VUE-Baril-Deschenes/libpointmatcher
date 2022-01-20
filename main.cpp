#include <iostream>
#include <fstream>
#include "pointmatcher/PointMatcher.h"
#include "pointmatcher/DataPointsFilters/SurfaceNormalGICPCovariance.h"

typedef PointMatcher<float> PM;

static double gicp_epsilon = 1e-3;

int main(int argc, char** argv)
{
    std::ifstream configFile("/home/dominic/sewervue_repos/libpointmatcher/test_config/realtime_post_filters.yaml");
//    PM::DataPoints reading = PM::DataPoints::load("/home/dominic/Desktop/1617238752456772327.vtk");
    PM::DataPoints reference = PM::DataPoints::load("/home/dominic/Desktop/1617238750956836700.vtk");
    PM::DataPointsFilters f(configFile);
    configFile.close();

    f.apply(reference);

    std::cout << reference.descriptorLabels << std::endl;
    std::cout << reference.getDescriptorViewByName("gicpCovariance").col(0) << std::endl;

    return 0;
}
