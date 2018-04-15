/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/FrictionDEM

================================================================================================= */

#ifndef FRICTIONFEM_GEOMETRY_CPP
#define FRICTIONFEM_GEOMETRY_CPP

// -------------------------------------------------------------------------------------------------

#include "Geometry.h"

// -------------------------------------------------------------------------------------------------

namespace FrictionDEM {

// -------------------------------------------------------------------------------------------------

inline Geometry::Geometry(const ColD &m, const MatD &x, const MatS &dofs) :
  m_x(x), m_m(m), m_dofs(dofs)
{
  // check input
  assert( m_x.size() == m_m.size()    );
  assert( m_x.rows() == m_dofs.rows() );
  assert( m_x.cols() == m_dofs.cols() );

  // allocate particle vectors
  m_v = MatD(m_x.rows(), m_x.cols());
  m_a = MatD(m_x.rows(), m_x.cols());

  // zero-initialize particle vectors
  m_v.setZero();
  m_a.setZero();
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::setSpring(const Spring &mat)
{
  m_spring = mat;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::x()
{
  return m_x;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::v()
{
  return m_v;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::a()
{
  return m_a;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::m()
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::f()
{
  // allocate force vector per particle
  MatD f(m_x.rows(), m_x.cols());

  // zero-initialize
  f.setZero();

  // evaluate constitutive models
  f += m_spring.force(m_x);

  return f;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_x()
{
}

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
