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
#include "IntensityNormalization.h"

// IntensityNormalizationDataPointsFilter
// Compute
template<typename T>
typename PointMatcher<T>::DataPoints 
IntensityNormalizationDataPointsFilter<T>::filter(const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void IntensityNormalizationDataPointsFilter<T>::inPlaceFilter(DataPoints& cloud)
{

	if (!cloud.descriptorExists("normals"))
		throw InvalidField("IntensityNormalizationDataPointsFilter: Error, cannot find normals in descriptors.");
	if (!cloud.descriptorExists("observationDirections"))
		throw InvalidField("IntensityNormalizationDataPointsFilter: Error, cannot find observation directions in descriptors.");

	cloud.allocateDescriptor("normalizedIntensity", 1);
	BOOST_AUTO(normalizedIntensities, cloud.getDescriptorViewByName("normalizedIntensity"));

	const BOOST_AUTO(normals, cloud.getDescriptorViewByName("normals"));
	const BOOST_AUTO(intensities, cloud.getDescriptorViewByName("intensity"));
	const BOOST_AUTO(observationDirections, cloud.getDescriptorViewByName("observationDirections"));
	assert(normals.rows() == observationDirections.rows());
	
	const unsigned int nbPts(cloud.getNbPoints());

    typename PointMatcher<T>::Matrix cosAngles (1, nbPts);

	for (unsigned int i = 0; i < nbPts; ++i)
	{
		// Check normal orientation
		const Vector vecP = observationDirections.col(i).normalized();
		const Vector vecN = normals.col(i);

        cosAngles(0,i) = vecP.dot(vecN);

        if (vecN.norm() < 0.9 || fabs(cosAngles(0,i)) < 0.1)
        {
            normalizedIntensities(0,i) = -1;
        }

        else
        {
            normalizedIntensities(0,i) = intensities(0,i) / cosAngles(0,i);
        }
	}
    // TODO: Re-code without the loop
}

template struct IntensityNormalizationDataPointsFilter<float>;
template struct IntensityNormalizationDataPointsFilter<double>;
