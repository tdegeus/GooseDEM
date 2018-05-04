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

// ======================= Trampoline for parent class - GooseDEM/Geometry.h =======================

class PyGeometry : public M::Geometry
{
public:
  // inherit the constructors
  using M::Geometry::Geometry;

  // trampoline
  ColD solve()                    override { PYBIND11_OVERLOAD_PURE( ColD, M::Geometry, solve); }
  void reset()                    override { PYBIND11_OVERLOAD_PURE( void, M::Geometry, reset); }
  bool stop(double tol)           override { PYBIND11_OVERLOAD_PURE( bool, M::Geometry, stop, tol); }
  void timestep(double dt)        override { PYBIND11_OVERLOAD_PURE( void, M::Geometry, timestep, dt); }
  MatD x() const                  override { PYBIND11_OVERLOAD_PURE( MatD, M::Geometry, x); }
  MatD v() const                  override { PYBIND11_OVERLOAD_PURE( MatD, M::Geometry, v); }
  MatD a() const                  override { PYBIND11_OVERLOAD_PURE( MatD, M::Geometry, a); }
  ColD dofs_v() const             override { PYBIND11_OVERLOAD_PURE( ColD, M::Geometry, dofs_v); }
  ColD dofs_a() const             override { PYBIND11_OVERLOAD_PURE( ColD, M::Geometry, dofs_a); }
  void set_x(const MatD &pvector) override { PYBIND11_OVERLOAD_PURE( void, M::Geometry, set_x, pvector); }
  void set_v(const ColD &dofval)  override { PYBIND11_OVERLOAD_PURE( void, M::Geometry, set_v, dofval); }
  void set_a(const ColD &dofval)  override { PYBIND11_OVERLOAD_PURE( void, M::Geometry, set_a, dofval); }
};

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

// ================================ GooseDEM - GooseDEM/Geometry.h =================================

py::class_<GooseDEM::Geometry, PyGeometry> geometry(m, "Geometry");

geometry
  // -
  .def(py::init<>())
  // -
  .def("solve"   , &M::Geometry::solve)
  .def("reset"   , &M::Geometry::reset)
  .def("stop"    , &M::Geometry::stop)
  .def("timestep", &M::Geometry::timestep)
  .def("x"       , &M::Geometry::x)
  .def("v"       , &M::Geometry::v)
  .def("a"       , &M::Geometry::a)
  .def("dofs_v"  , &M::Geometry::dofs_v)
  .def("dofs_a"  , &M::Geometry::dofs_a)
  .def("set_x"   , &M::Geometry::set_x)
  .def("set_v"   , &M::Geometry::set_v)
  .def("set_v"   , &M::Geometry::set_v)
  // print to screen
  .def("__repr__",
    [](const GooseDEM::Geometry &a){ return "<GooseDEM.Geometry>"; }
  );

// ============================ GooseDEM - GooseDEM/GeometryFriction.h =============================

py::class_<GooseDEM::GeometryFriction>(m, "GeometryFriction", geometry)
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
  .def("fix_v", &M::GeometryFriction::fix_v)
  .def("fext", &M::GeometryFriction::fext)
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
  .def("set_v", py::overload_cast<cColD &>(&M::GeometryFriction::set_v))
  .def("set_a", py::overload_cast<cColD &>(&M::GeometryFriction::set_a))
  .def("set_x",                            &M::GeometryFriction::set_x )
  .def("set_v", py::overload_cast<cMatD &>(&M::GeometryFriction::set_v))
  .def("set_a", py::overload_cast<cMatD &>(&M::GeometryFriction::set_a))
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

