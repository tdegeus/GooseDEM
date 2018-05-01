/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_GEOMETRY_FRICTION_H
#define GOOSEDEM_GEOMETRY_FRICTION_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class GeometryFriction : public Geometry
{
private:

  // conversion vector
  Vector m_vec;

  // particles
  MatD m_fext;   // external force [N, ndim]
  MatD m_x;      // position       [N, ndim]
  MatD m_v;      // velocity       [N, ndim]
  MatD m_a;      // acceleration   [N, ndim]
  ColD m_m;      // mass           [N, ndim]
  MatS m_dofs;   // DOF-number     [N, ndim]

  // prescribed DOFs
  ColS m_iip;    // DOF-numbers         [np]
  ColD m_vp;     // prescribed velocity [np]

  // DOF values
  ColD m_M;      // mass            [ndof]
  ColD m_Minv;   // inverse of mass [ndof]

  // time & convergence check
  double   m_t;
  StopList m_stop;

  // dimensions
  size_t m_N;    // number of particles
  size_t m_ndim; // number of spatial dimensions
  size_t m_ndof; // number of DOFs

  // constitutive models
  Spring  m_spring;
  Dashpot m_dashpot;

public:

  // constructor
  GeometryFriction(const ColD &m, const MatD &x, const MatS &dofs);

  // append constitutive models
  void set(const Spring  &mat);
  void set(const Dashpot &mat);

  // set fixed velocity
  void fix_v(const ColS &iip, const ColD &vp);

  // set external force
  void fext(const MatD &pvector);

  // solve for DOF-accelerations [ndof]
  ColD solve() override;

  // reset residuals, check for convergence
  void reset() override;
  bool stop(double tol) override;

  // process time-step
  void timestep(double dt) override;

  // return particle vectors [N, ndim]
  MatD x() const override;
  MatD v() const override;
  MatD a() const override;
  MatD f() const;
  ColD m() const;

  // return DOF values [ndof]
  ColD dofs_v() const override;
  ColD dofs_a() const override;
  ColD dofs_f() const;
  ColD dofs_m() const;

  // overwrite particle vectors [N, ndim]
  void set_x(const MatD &pvector) override;
  void set_v(const MatD &pvector);
  void set_a(const MatD &pvector);

  // overwrite particle vectors, reconstructed from DOF values
  void set_v(const ColD &dofval) override; // == set_v(asParticle(V))
  void set_a(const ColD &dofval) override; // == set_a(asParticle(A))

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
