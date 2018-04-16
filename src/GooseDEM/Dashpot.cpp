/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_DASHPOT_CPP
#define GOOSEDEM_DASHPOT_CPP

// -------------------------------------------------------------------------------------------------

#include "Dashpot.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline Dashpot::Dashpot(const MatS &particles, const ColD &eta) :
  m_particles(particles), m_eta(eta)
{
  // check input
  assert( m_particles.rows() == m_eta.size() );
}

// -------------------------------------------------------------------------------------------------

inline MatD Dashpot::force(const MatD &V) const
{
  // dimensions
  auto n    = V.rows(); // number of particles
  auto ndim = V.cols(); // number of dimensions

  // force per particle
  // - allocate
  MatD F(n, ndim);
  // - zero-initialize
  F.setZero();

  // local variables
  cppmat::cartesian::vector<double> vi(ndim); // velocity of particle "i"
  cppmat::cartesian::vector<double> vj(ndim); // velocity of particle "j"
  cppmat::cartesian::vector<double> dv(ndim); // velocity difference
  cppmat::cartesian::vector<double> f (ndim); // force vector

  // loop over all dashpots
  for ( auto p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    auto i = m_particles(p,0);
    auto j = m_particles(p,1);
    // - copy the particles' velocities to the vectors "vi" and "vj"
    std::copy(V.data()+i*ndim, V.data()+(i+1)*ndim, vi.data());
    std::copy(V.data()+j*ndim, V.data()+(j+1)*ndim, vj.data());
    // - compute the velocity difference vector
    dv = vj - vi;
    // - compute the force vector
    f = m_eta(p) * dv;
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

inline ColS Dashpot::coordination(const MatD &X) const
{
  // coordination per particle
  // - allocate
  ColS C(X.rows());
  // - zero-initialize
  C.setZero();

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

}

// -------------------------------------------------------------------------------------------------

#endif
