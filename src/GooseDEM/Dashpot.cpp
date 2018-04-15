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
  // number of spatial dimensions
  int nd = V.cols();

  // force per particle
  // - allocate
  MatD F(V.rows(), nd);
  // - zero-initialize
  F.setZero();

  // local variables
  cppmat::cartesian::vector<double> vi(nd); // velocity of particle "i"
  cppmat::cartesian::vector<double> vj(nd); // velocity of particle "j"
  cppmat::cartesian::vector<double> dv(nd); // velocity difference
  cppmat::cartesian::vector<double> f (nd); // force vector

  // loop over all dashpots
  for ( size_t p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    size_t i = m_particles(p,0);
    size_t j = m_particles(p,1);
    // - copy the particles' velocities to the vectors "vi" and "vj"
    std::copy(V.data()+i*nd, V.data()+(i+1)*nd, vi.data());
    std::copy(V.data()+j*nd, V.data()+(j+1)*nd, vj.data());
    // - compute the velocity difference vector
    dv = vj - vi;
    // - compute the force vector
    f = m_eta(p) * dv;
    // - assemble the force to the particles
    for ( size_t d = 0 ; d < nd ; ++d )
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
  for ( size_t p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    size_t i = m_particles(p,0);
    size_t j = m_particles(p,1);
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
