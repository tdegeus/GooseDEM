/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_SPRING_CPP
#define GOOSEDEM_SPRING_CPP

// -------------------------------------------------------------------------------------------------

#include "Spring.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline Spring::Spring(const MatS &particles, const ColD &k, const ColD &D0) :
  m_particles(particles), m_k(k), m_D0(D0)
{
  // check input
  assert( m_particles.rows() == m_k .size() );
  assert( m_particles.rows() == m_D0.size() );
}

// -------------------------------------------------------------------------------------------------

inline MatD Spring::force(const MatD &X) const
{
  // number of spatial dimensions
  int nd = X.cols();

  // force per particle
  // - allocate
  MatD F(X.rows(), nd);
  // - zero-initialize
  F.setZero();

  // local variables
  cppmat::cartesian::vector<double> xi(nd); // position of particle "i"
  cppmat::cartesian::vector<double> xj(nd); // position of particle "j"
  cppmat::cartesian::vector<double> dx(nd); // position difference
  cppmat::cartesian::vector<double> f (nd); // force vector
  double D; // distance in 'local coordinates'

  // loop over all springs
  for ( size_t p = 0 ; p < m_particles.rows() ; ++p )
  {
    // - extract particle numbers
    size_t i = m_particles(p,0);
    size_t j = m_particles(p,1);
    // - copy the particles' positions to the vectors "xi" and "xj"
    std::copy(X.data()+i*nd, X.data()+(i+1)*nd, xi.data());
    std::copy(X.data()+j*nd, X.data()+(j+1)*nd, xj.data());
    // - compute the position difference vector
    dx = xj - xi;
    // - compute the current length
    D = dx.length();
    // - compute the force vector, by comparing to the spring's relaxed length
    f = m_k(p) * (D - m_D0(p)) * dx/D;
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

inline ColS Spring::coordination(const MatD &X) const
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
