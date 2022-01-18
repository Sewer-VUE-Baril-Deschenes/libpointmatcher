#ifndef LIBPOINTMATCHER_GICP_H
#define LIBPOINTMATCHER_GICP_H


#include <PointMatcher.h>

template<typename T>
struct GICPErrorMinimizer: PointMatcher<T>::ErrorMinimizer
{
	typedef PointMatcherSupport::Parametrizable Parametrizable;
	typedef PointMatcherSupport::Parametrizable P;
	typedef Parametrizable::Parameters Parameters;
	typedef Parametrizable::ParametersDoc ParametersDoc;

	typedef typename PointMatcher<T>::ErrorMinimizer::ErrorElements ErrorElements;
	typedef typename PointMatcher<T>::TransformationParameters TransformationParameters;

	virtual inline const std::string name()
	{
		return "GICPErrorMinimizer";
	}

	inline static const std::string description()
	{
		return "Generalized ICP error. Per \\cite{Segal2009}.";
	}

	inline static const ParametersDoc availableParameters()
	{
		return {};
	}

	GICPErrorMinimizer(const Parameters& params = Parameters());
	GICPErrorMinimizer(const ParametersDoc paramsDoc, const Parameters& params);
	virtual TransformationParameters compute(const ErrorElements& mPts);
};


#endif
