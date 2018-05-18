/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "Friction.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = GooseDEM;
namespace E  = GooseDEM::Ext::Friction;

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

// =================================================================================================

PYBIND11_MODULE(GooseDEM_Ext_Friction, m) {

m.doc() = "Friction simulation";

// =========================== GooseDEM/Ext/Friction/PotentialAdhsion.h ============================

py::class_<E::PotentialAdhesion>(m, "PotentialAdhesion")
  // constructor
  .def(
    py::init<cMatS &, cColD &, cColD &, cColD &, cColD &>(),
    "PotentialAdhesion",
    py::arg("particles"),
    py::arg("k"),
    py::arg("b"),
    py::arg("r0"),
    py::arg("e")
  )
  // methods
  .def("force"       , &E::PotentialAdhesion::force       )
  .def("coordination", &E::PotentialAdhesion::coordination)
  .def("potential"   , &E::PotentialAdhesion::potential   )
  // print to screen
  .def("__repr__",
    [](const E::PotentialAdhesion &a){ return "<GooseDEM_Ext_Friction.PotentialAdhesion>"; }
  );

// =============================== GooseDEM/Ext/Friction/Geometry.h ================================

py::class_<E::Geometry, M::Geometry>(m, "Geometry")
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
  .def("set", py::overload_cast<const M::Spring            &>(&E::Geometry::set))
  .def("set", py::overload_cast<const M::Dashpot           &>(&E::Geometry::set))
  .def("set", py::overload_cast<const E::PotentialAdhesion &>(&E::Geometry::set))
  // -
  .def("fix_v"    , &E::Geometry::fix_v   )
  .def("set_fext" , &E::Geometry::set_fext)
  .def("fext"     , &E::Geometry::fext    )
  .def("fres"     , &E::Geometry::fres    )
  // -
  .def("x"           , &E::Geometry::x           )
  .def("v"           , &E::Geometry::v           )
  .def("a"           , &E::Geometry::a           )
  .def("f"           , &E::Geometry::f           )
  .def("m"           , &E::Geometry::m           )
  .def("coordination", &E::Geometry::coordination)
  // -
  .def("dofs_v", &E::Geometry::dofs_v)
  .def("dofs_a", &E::Geometry::dofs_a)
  .def("dofs_f", &E::Geometry::dofs_f)
  .def("dofs_m", &E::Geometry::dofs_m)
  // -
  .def("set_v", py::overload_cast<cColD &>(&E::Geometry::set_v))
  .def("set_a", py::overload_cast<cColD &>(&E::Geometry::set_a))
  .def("set_x",                            &E::Geometry::set_x )
  .def("set_v", py::overload_cast<cMatD &>(&E::Geometry::set_v))
  .def("set_a", py::overload_cast<cMatD &>(&E::Geometry::set_a))
  // print to screen
  .def("__repr__",
    [](const E::Geometry &a){ return "<GooseDEM_Ext_Friction.Geometry>"; }
  );

// =================================================================================================

}
