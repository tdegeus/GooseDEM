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
typedef GooseDEM::ColD ColD;
typedef GooseDEM::ColS ColS;
typedef GooseDEM::MatD MatD;
typedef GooseDEM::MatS MatS;

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

// ================================= GooseDEM - GooseDEM/Geometry.h ==================================

py::class_<GooseDEM::Geometry>(m, "Geometry")
  // constructor
  .def(
    py::init<cColD &, cMatD &, cMatS &>(),
    "Geometry",
    py::arg("m"),
    py::arg("x"),
    py::arg("dofs")
  )
  // methods
  // -
  .def("set", py::overload_cast<const M::Spring  &>(&M::Geometry::set))
  .def("set", py::overload_cast<const M::Dashpot &>(&M::Geometry::set))
  // -
  .def("x", &M::Geometry::x)
  .def("v", &M::Geometry::v)
  .def("a", &M::Geometry::a)
  .def("f", &M::Geometry::f)
  .def("m", &M::Geometry::m)
  // -
  .def("dofs_x", &M::Geometry::dofs_x)
  .def("dofs_v", &M::Geometry::dofs_v)
  .def("dofs_a", &M::Geometry::dofs_a)
  .def("dofs_f", &M::Geometry::dofs_f)
  .def("dofs_m", &M::Geometry::dofs_m)
  // -
  .def("set_x", &M::Geometry::set_x)
  .def("set_v", &M::Geometry::set_v)
  .def("set_a", &M::Geometry::set_a)
  // -
  .def("asDofs"      , &M::Geometry::asDofs      )
  .def("assembleDofs", &M::Geometry::assembleDofs)
  .def("asParticle"  , &M::Geometry::asParticle  )
  // print to screen
  .def("__repr__",
    [](const GooseDEM::Geometry &a){ return "<GooseDEM.Geometry>"; }
  );

// -------------------------------------------------------------------------------------------------

// m.def("dump",&M::dump,
//   "Dump to text file",
//   py::arg("fname"),
//   py::arg("matrix")
// );

// =================================== GooseDEM - GooseDEM/Sim.h ===================================

m.def("velocityVerlet",&M::velocityVerlet,
  "evaluate one time step",
  py::arg("geometry"),
  py::arg("dt"),
  py::arg("iip"),
  py::arg("v_p")
);

// -------------------------------------------------------------------------------------------------

m.def("quasiStaticVelocityVerlet",&M::quasiStaticVelocityVerlet,
  "iterate until all particles have come to a rest",
  py::arg("geometry"),
  py::arg("dt"),
  py::arg("iip"),
  py::arg("v_p"),
  py::arg("norm"),
  py::arg("ncheck")=20
);

// =================================================================================================

}

