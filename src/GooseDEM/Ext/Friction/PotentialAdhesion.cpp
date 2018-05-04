/* =================================================================================================

Contributors:

  Tianxia Ma
  Wencheng Ji

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_EXT_FRICTION_POTENTIALADHESION_CPP
#define GOOSEDEM_EXT_FRICTION_POTENTIALADHESION_CPP

// -------------------------------------------------------------------------------------------------

#include "PotentialAdhesion.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {
namespace Ext {
namespace Friction {

// -------------------------------------------------------------------------------------------------

inline PotentialAdhesion::PotentialAdhesion(const MatS &particles, const ColD &k, const ColD &D0) :
  m_particles(particles), m_k(k), m_D0(D0)
{
  // check input
  assert( m_particles.rows() == m_k .size() );
  assert( m_particles.rows() == m_D0.size() );
}

// -------------------------------------------------------------------------------------------------

inline MatD PotentialAdhesion::force(const MatD &X) const
{
  // dimensions
  auto n    = X.rows(); // number of particles
  auto ndim = X.cols(); // number of dimensions

  // zero-initialize force per particle
  MatD F = MatD::Zero(n, ndim);

  // local variables
  cppmat::cartesian::vector<double> xi(ndim); // position of particle "i"
  cppmat::cartesian::vector<double> xj(ndim); // position of particle "j"
  cppmat::cartesian::vector<double> dx(ndim); // position difference
  cppmat::cartesian::vector<double> f (ndim); // force vector
  double D; // distance in 'local coordinates'

  // loop over all springs
  for ( auto p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    auto i = m_particles(p,0);
    auto j = m_particles(p,1);
    // - copy the particles' positions to the vectors "xi" and "xj"
    std::copy(X.data()+i*ndim, X.data()+(i+1)*ndim, xi.data());
    std::copy(X.data()+j*ndim, X.data()+(j+1)*ndim, xj.data());
    // - compute the position difference vector
    dx = xj - xi;
    // - compute the current length
    D = dx.length();
    // - compute the force vector, by comparing to the spring's relaxed length
    f = m_k(p) * (D - m_D0(p)) * dx/D;
    // - assemble the force to the particles
    for ( auto d = 0 ; d < ndim ; ++d )
    {
      F(i,d) += f(d);
      F(j,d) -= f(d);
    }
  }

  return F;
}

// -------------------------------------------------------------------------------------------------

inline ColS PotentialAdhesion::coordination(const MatD &X) const
{
  // zero-initialize coordination per particle
  ColS C = ColS::Zero(X.rows());

  // loop over all springs
  for ( auto p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    auto i = m_particles(p,0);
    auto j = m_particles(p,1);
    // - add to coordination
    C(i) += 1;
    C(j) += 1;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
