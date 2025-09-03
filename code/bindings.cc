#include <pybind11/pybind11.h>
#include "CohesiveCrack.cc"

namespace py = pybind11;

PYBIND11_MODULE(CohesiveCrack, m) {
    m.doc() = "Cohesive crack stress field analysis";

    m.def("delta_sigma_xy", &StressAnalysis::delta_sigma_xy,
          "Compute shear stress component",
          py::arg("x"), py::arg("y"), py::arg("X_c"),
          py::arg("C_f"), py::arg("C_s"), py::arg("C_d"),
          py::arg("nu"), py::arg("Gamma"), py::arg("E"));

    m.def("delta_sigma_xx", &StressAnalysis::delta_sigma_xx,
          "Compute normal stress component",
          py::arg("x"), py::arg("y"), py::arg("X_c"),
          py::arg("C_f"), py::arg("C_s"), py::arg("C_d"),
          py::arg("nu"), py::arg("Gamma"), py::arg("E"));
}