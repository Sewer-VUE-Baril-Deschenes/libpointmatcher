#include "gicp.h"
#include "pointmatcher/ErrorMinimizersImpl.h"

namespace python
{
	namespace errorminimizers
	{
		void pybindGICP(py::module& p_module)
		{
			using GICPErrorMinimizer = ErrorMinimizersImpl<ScalarType>::GICPErrorMinimizer;
			py::class_<GICPErrorMinimizer, std::shared_ptr<GICPErrorMinimizer>, ErrorMinimizer>(p_module, "GICPErrorMinimizer")
					.def(py::init<const Parameters&>(), py::arg("params") = Parameters())
					.def(py::init<const ParametersDoc, const Parameters&>(), py::arg("paramsDoc"), py::arg("params") = Parameters())

					.def_static("description", &GICPErrorMinimizer::description)
					.def_static("availableParameters", &GICPErrorMinimizer::availableParameters)

					.def("name", &GICPErrorMinimizer::name)
					.def("compute", &GICPErrorMinimizer::compute, py::arg("mPts"));
		}
	}
}
