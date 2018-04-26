/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "../src/GooseDEM/GooseDEM.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = GooseDEM;

// abbreviate types(s)
// - Eigen types
typedef GooseDEM::ColD ColD;
typedef GooseDEM::ColS ColS;
typedef GooseDEM::MatD MatD;
typedef GooseDEM::MatS MatS;
// - const arrays
typedef const GooseDEM::ColD cColD;
typedef const GooseDEM::ColS cColS;
typedef const GooseDEM::MatD cMatD;
typedef const GooseDEM::MatS cMatS;

// =========================================== GooseDEM ============================================

PYBIND11_MODULE(GooseDEM, m) {

m.doc() = "Simple DEM simulation";

// ================================= GooseDEM - GooseDEM/Spring.h ==================================

py::class_<GooseDEM::Spring>(m, "Spring")
  // constructor
  .def(
    py::init<cMatS &, cColD &, cColD &>(),
    "Spring",
    py::arg("particles"),
    py::arg("k"),
    py::arg("D0")
  )
  // methods
  .def("force"       , &M::Spring::force       )
  .def("coordination", &M::Spring::coordination)
  // print to screen
  .def("__repr__",
    [](const GooseDEM::Spring &a){ return "<GooseDEM.Spring>"; }
  );

// ================================= GooseDEM - GooseDEM/Dashpot.h ==================================

py::class_<GooseDEM::Dashpot>(m, "Dashpot")
  // constructor
  .def(
    py::init<cMatS &, cColD &>(),
    "Dashpot",
    py::arg("particles"),
    py::arg("eta")
  )
  // methods
  .def("force"       , &M::Dashpot::force       )
  .def("coordination", &M::Dashpot::coordination)
  // print to screen
  .def("__repr__",
    [](const GooseDEM::Dashpot &a){ return "<GooseDEM.Dashpot>"; }
  );

// ============================ GooseDEM - GooseDEM/GeometryFriction.h =============================

py::class_<GooseDEM::GeometryFriction>(m, "GeometryFriction")
  // constructor
  .def(
    py::init<cColD &, cMatD &, cMatS &>(),
    "GeometryFriction",
    py::arg("m"),
    py::arg("x"),
    py::arg("dofs")
  )
  // methods
  // -
  .def("set", py::overload_cast<const M::Spring  &>(&M::GeometryFriction::set))
  .def("set", py::overload_cast<const M::Dashpot &>(&M::GeometryFriction::set))
  // -
  .def("x", &M::GeometryFriction::x)
  .def("v", &M::GeometryFriction::v)
  .def("a", &M::GeometryFriction::a)
  .def("f", &M::GeometryFriction::f)
  .def("m", &M::GeometryFriction::m)
  // -
  .def("dofs_v", &M::GeometryFriction::dofs_v)
  .def("dofs_a", &M::GeometryFriction::dofs_a)
  .def("dofs_f", &M::GeometryFriction::dofs_f)
  .def("dofs_m", &M::GeometryFriction::dofs_m)
  // -
  .def("set_v"     , &M::GeometryFriction::set_v)
  .def("set_a"     , &M::GeometryFriction::set_a)
  .def("set_x"     , &M::GeometryFriction::set_x )
  .def("set_v_dofs", &M::GeometryFriction::set_v_dofs)
  .def("set_a_dofs", &M::GeometryFriction::set_a_dofs)
  // print to screen
  .def("__repr__",
    [](const GooseDEM::GeometryFriction &a){ return "<GooseDEM.GeometryFriction>"; }
  );

// ================================== GooseDEM - GooseDEM/Write.h ==================================

m.def("dump", py::overload_cast<const std::string &, cColD &>(&M::dump),
  "Dump to text file",
  py::arg("fname"),
  py::arg("matrix")
);

// -------------------------------------------------------------------------------------------------

m.def("dump", py::overload_cast<const std::string &, cMatD &>(&M::dump),
  "Dump to text file",
  py::arg("fname"),
  py::arg("matrix")
);

// ============================= GooseDEM - GooseDEM/TimeIntegration.h =============================

m.def("velocityVerlet", &M::velocityVerlet,
  "evaluate one time step",
  py::arg("geometry"),
  py::arg("dt")
);

// -------------------------------------------------------------------------------------------------

m.def("quasiStaticVelocityVerlet", &M::quasiStaticVelocityVerlet,
  "iterate until all particles have come to a rest",
  py::arg("geometry"),
  py::arg("dt"),
  py::arg("tol")
);

// =================================================================================================

}

