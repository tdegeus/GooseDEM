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
#include <math.h>

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {
namespace Ext {
namespace Friction {

// -------------------------------------------------------------------------------------------------

inline PotentialAdhesion::PotentialAdhesion(const MatS &particles, const ColD &k, const ColD &b, const ColD &r0, const ColD &e) :
  m_particles(particles), m_k(k), m_b(b), m_r0(r0), m_e(e)
{
  // check input
  assert( m_particles.rows() == m_k .size() );
  assert( m_particles.rows() == m_b .size() );
  assert( m_particles.rows() == m_r0.size() );
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

  // loop over all interacted particle pairs
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
    // - compute the force vector, by comparing to the interacted particles' initial length
    if ( D <= m_r0(p) ) {
        f = 12 / D * ( std::pow( m_r0(p) / D, 6 ) - std::pow( m_r0(p) / D, 12 ) ); // - Lennard-Jones force
    }
    else {
        f = (std::pow( m_k(p), 2 ) * m_b(p) * (std::pow( D - m_r0(p), 2) ) + 2 * m_k(p) * (D - m_r0(p))) * exp(- m_k(p) * m_b(p) * (D - m_r0(p))) * dx/D; // - Innovative force
    }

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

  // loop over all interacted particle pairs
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

inline ColD PotentialAdhesion::potential(const MatD &X) const
{
  // dimensions
  auto n    = X.rows(); // number of particles
  auto ndim = X.cols(); // number of dimensions

  // zero-initialize force per particle
  MatD V = MatD::Zero(n, 1);

  // local variables
  cppmat::cartesian::vector<double> xi(ndim); // position of particle "i"
  cppmat::cartesian::vector<double> xj(ndim); // position of particle "j"
  cppmat::cartesian::vector<double> dx(ndim); // position difference
  // cppmat::cartesian::vector<double> v (ndim); // potential vector
  double D; // distance in 'local coordinates'

  // loop over all interacted particle pairs
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
    // - compute the force vector, by comparing to the interacted particles' initial length
    if ( D <= m_r0(p) ) {
        V(p, 0) = m_e(p) * ( std::pow( m_r0(p) / D, 12 ) - 2 * std::pow( m_r0(p) / D, 6 ) );
    }
    else {
        V(p, 0) = - m_k(p) * std::pow( D - m_r0(p) + 2 / (m_k(p) * m_b(p)), 2 ) * exp( - m_k(p) * m_b(p) * (D - m_r0(p)) ) + 4 / (m_k(p) * pow(m_b(p), 2)) - m_e(p);
    }

    // // - assemble the force to the particles
    // for ( auto d = 0 ; d < ndim ; ++d )
    // {
    //   V(i,d) += v(d);
    //   V(j,d) -= v(d);
    // }
  }

  return V;
}

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
