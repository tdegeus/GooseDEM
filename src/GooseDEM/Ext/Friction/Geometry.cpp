/* =================================================================================================

Contributors:

  Tianxia Ma
  Wencheng Ji

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_EXT_FRICITION_GEOMETY_CPP
#define GOOSEDEM_EXT_FRICITION_GEOMETY_CPP

// -------------------------------------------------------------------------------------------------

#include "Geometry.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {
namespace Ext {
namespace Friction {

// -------------------------------------------------------------------------------------------------

inline Geometry::Geometry(const ColD &m, const MatD &x, const MatS &dofs) :
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

inline void Geometry::set(const Spring &mat)
{
  m_spring = mat;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set(const Dashpot &mat)
{
  m_dashpot = mat;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set(const PotentialAdhesion &mat)
{
  m_potentialadhesion = mat;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::fix_v(const ColS &iip, const ColD &vp)
{
  m_iip = iip;
  m_vp  = vp;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::fext(const MatD &pvector)
{
  assert( pvector.rows() == m_x.rows() );
  assert( pvector.cols() == m_x.cols() );

  m_fext = pvector;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::solve()
{
  return m_Minv.cwiseProduct( m_vec.assembleDofs( f() + m_fext ) );
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::reset()
{
  m_stop.reset();
}

// -------------------------------------------------------------------------------------------------

inline bool Geometry::stop(double tol)
{
  // external and internal force
  ColD Fext = m_vec.assembleDofs(m_fext);
  ColD Fint = m_vec.assembleDofs(f());

  // sum of absolute
  double res  = Fint.cwiseAbs().sum();
  double fext = Fext.cwiseAbs().sum();

  // normalize
  if ( fext != 0 ) res /= fext;

  // evaluate stopping criterion
  return m_stop.stop(res, tol);
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::timestep(double dt)
{
  m_t += dt;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::x() const
{
  return m_x;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::v() const
{
  return m_v;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::a() const
{
  return m_a;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::m() const
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::f() const
{
  // zero-initialize force vector per particle
  MatD f = MatD::Zero(m_N, m_ndim);

  // evaluate constitutive models
  f += m_spring           .force(m_x);
  f += m_dashpot          .force(m_v);
  f += m_potentialadhesion.force(m_v);

  return f;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_v() const
{
  return m_vec.asDofs(m_v);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_a() const
{
  return m_vec.asDofs(m_a);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_m() const
{
  return m_vec.asDofs(m_m);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_f() const
{
  return m_vec.assembleDofs(f());
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_x(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_x = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_v(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_v = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_a(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_a = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_v(const ColD &Vin)
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

inline void Geometry::set_a(const ColD &Ain)
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

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
