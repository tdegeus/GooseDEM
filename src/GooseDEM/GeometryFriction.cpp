/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_GEOMETRY_FRICTION_CPP
#define GOOSEDEM_GEOMETRY_FRICTION_CPP

// -------------------------------------------------------------------------------------------------

#include "GeometryFriction.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline GeometryFriction::GeometryFriction(const ColD &m, const MatD &x, const MatS &dofs) :
  m_x(x), m_m(m), m_dofs(dofs)
{
  // extract dimensions
  m_N    = static_cast<size_t>(m_dofs.rows());
  m_ndim = static_cast<size_t>(m_dofs.cols());
  m_ndof = static_cast<size_t>(m_dofs.maxCoeff() + 1);

  // check input
  assert( m_x.rows() == m_N    );
  assert( m_x.cols() == m_ndim );
  assert( m_m.size() == m_N    );

  // conversion vector
  m_vec  = Vector(dofs);

  // zero-initialize particle vectors
  m_v    = MatD::Zero(m_x.rows(), m_x.cols());
  m_a    = MatD::Zero(m_x.rows(), m_x.cols());

  // zero-initialize boundary conditions
  m_iip  = ColS();
  m_vp   = ColD();
  m_fext = MatD::Zero(m_x.rows(), m_x.cols());

  // zero-initialize time
  m_t    = 0.0;

  // compute (inverse of) DOF masses
  m_M    = m_vec.asDofs(m_m);
  m_Minv = m_M.cwiseInverse();
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set(const Spring &mat)
{
  m_spring = mat;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set(const Dashpot &mat)
{
  m_dashpot = mat;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::fix_v(const ColS &iip, const ColD &vp)
{
  m_iip = iip;
  m_vp  = vp;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::fext(const MatD &pvector)
{
  assert( pvector.rows() == m_x.rows() );
  assert( pvector.cols() == m_x.cols() );

  m_fext = pvector;
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::solve()
{
  return m_Minv.cwiseProduct( m_vec.assembleDofs( f() + m_fext ) );
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::reset()
{
  m_stop.reset();
}

// -------------------------------------------------------------------------------------------------

inline bool GeometryFriction::stop(double tol)
{
  // parameters
  double res = 0.0;
  ColD   V   = dofs_v();

  // exclude fixed velocities
  for ( auto i = 0 ; i < m_iip.size() ; ++i ) V(m_iip(i)) = 0.0;

  // compute kinetic energy
  for ( auto i = 0 ; i < V.size() ; ++i )
      res += 0.5 * m_M(i) * std::pow(V(i),2.0);

  return m_stop.stop(res, tol);
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::timestep(double dt)
{
  m_t += dt;
}

// -------------------------------------------------------------------------------------------------

inline MatD GeometryFriction::x() const
{
  return m_x;
}

// -------------------------------------------------------------------------------------------------

inline MatD GeometryFriction::v() const
{
  return m_v;
}

// -------------------------------------------------------------------------------------------------

inline MatD GeometryFriction::a() const
{
  return m_a;
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::m() const
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

inline MatD GeometryFriction::f() const
{
  // zero-initialize force vector per particle
  MatD f = MatD::Zero(m_N, m_ndim);

  // evaluate constitutive models
  f += m_spring .force(m_x);
  f += m_dashpot.force(m_x);

  return f;
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::dofs_v() const
{
  return m_vec.asDofs(m_v);
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::dofs_a() const
{
  return m_vec.asDofs(m_a);
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::dofs_m() const
{
  return m_vec.asDofs(m_m);
}

// -------------------------------------------------------------------------------------------------

inline ColD GeometryFriction::dofs_f() const
{
  return m_vec.assembleDofs(f());
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set_x(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_x = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set_v(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_v = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set_a(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_a = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set_v_dofs(const ColD &Vin)
{
  // check input
  assert( Vin.size() == m_ndof );

  // copy input
  ColD V = Vin;

  // apply boundary conditions
  for ( auto i = 0 ; i < m_iip.size() ; ++i ) V(m_iip(i)) = m_vp(i);

  // reconstruct and save nodal quantities
  m_v = m_vec.asParticle(V);
}

// -------------------------------------------------------------------------------------------------

inline void GeometryFriction::set_a_dofs(const ColD &Ain)
{
  // check input
  assert( Ain.size() == m_ndof );

  // copy input
  ColD A = Ain;

  // apply boundary conditions
  for ( auto i = 0 ; i < m_iip.size() ; ++i ) A(m_iip(i)) = 0.0;

  // reconstruct and save nodal quantities
  m_v = m_vec.asParticle(A);
}

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
