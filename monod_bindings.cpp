#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "monod_model.hpp"

namespace py = pybind11;

PYBIND11_MODULE(monod_cpp, m) {
    m.doc() = "Pybind11 interface for Monod bioreaction model";

    m.def("dX_dt", &monod::dX_dt, py::arg("biomass"), py::arg("substrate"),
          "Compute biomass derivative for Monod kinetics");
    m.def("dS_dt", &monod::dS_dt, py::arg("biomass"), py::arg("substrate"),
          "Compute substrate derivative for Monod kinetics");
    m.def(
        "integrate",
        [](double x0, double s0, double t_max, double dt) {
            const monod::SimulationResult result = monod::integrate(x0, s0, t_max, dt);
            return py::make_tuple(result.time_points, result.biomass, result.substrate);
        },
        py::arg("x0"), py::arg("s0"), py::arg("t_max"), py::arg("dt"),
        "Integrate Monod model with explicit Euler method");
}
