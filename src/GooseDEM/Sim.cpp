/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_SIM_CPP
#define GOOSEDEM_SIM_CPP

// -------------------------------------------------------------------------------------------------

#include "Sim.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline void velocityVerlet(Geometry &g, double dt, const ColS &iip, const ColD &v_p)
{
  // local variables and history
  ColD V;
  ColD V_n  = g.dofs_v();
  ColD A;
  ColD A_n  = g.dofs_a();
  ColD M    = g.dofs_m();
  ColD Minv = M.cwiseInverse();

  // (1) new positions: x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

  // (2a) estimate new velocities
  // - update velocities (DOFs)
  V = V_n + dt * A_n;
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - update velocities (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );

  // (2b) new velocities
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - update velocities (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );

  // (3) new accelerations
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - enforce prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) A(iip(i)) = 0.0;
  // - update accelerations (particles)
  g.set_a( g.asParticle(A) );
}

// -------------------------------------------------------------------------------------------------

inline size_t quasiStaticVelocityVerlet(
  Geometry &g, double dt, const ColS &iip, const ColD &v_p, double norm, size_t ncheck)
{
  // local variables
  ColD   V, V_n, A, A_n;
  ColD   M     = g.dofs_m();
  ColD   Minv  = M.cwiseInverse();
  size_t iiter = 0;
  bool   stop;
  double new_res;

  // list with residuals
  std::vector<double> res(ncheck, std::numeric_limits<double>::infinity());

  // 'iterate'
  while ( true )
  {
    // update history
    // --------------

    // update iteration counter counter
    ++iiter;

    // store history
    ColD V_n = g.dofs_v();
    ColD A_n = g.dofs_a();

    // velocity Verlet
    // ---------------

    // (1) new positions: x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
    g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

    // (2a) estimate new velocities
    // - update velocities (DOFs)
    V = V_n + dt * A_n;
    // - set prescribed velocities
    for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
    // - update velocities (particles)
    g.set_v( g.asParticle(V) );
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( - g.dofs_f() );
    // - update velocities (DOFs)
    V = V_n + .5 * dt * ( A_n + A );
    // - set prescribed velocities
    for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
    // - update velocities (particles)
    g.set_v( g.asParticle(V) );

    // (2b) new velocities
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( - g.dofs_f() );
    // - update velocities (DOFs)
    V = V_n + .5 * dt * ( A_n + A );
    // - set prescribed velocities
    for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = v_p(i);
    // - update velocities (particles)
    g.set_v( g.asParticle(V) );

    // (3) new accelerations
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( - g.dofs_f() );
    // - enforce prescribed velocities
    for ( int i = 0 ; i < iip.size() ; ++i ) A(iip(i)) = 0.0;
    // - update accelerations (particles)
    g.set_a( g.asParticle(A) );

    // check for convergence
    // ---------------------

    // compute total kinetic energy
    // - ignore prescribed DOFs
    for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = 0.0;
    // - zero-initialize
    new_res = 0.0;
    // - compute
    for ( size_t i = 0 ; i < static_cast<size_t>(V.size()) ; ++i )
      new_res += 0.5 * M(i) * std::pow(V(i),2.0);

    // move residual one place back (forgetting the first)
    for ( size_t i = 1 ; i < ncheck ; ++i )
      res[i-1] = res[i];

    // add new residual to the end
    res[ncheck-1] = new_res;

    // check for convergence: all residuals should be below the norm, and shouldn't be increasing
    // - initialize
    stop = true;
    // - evaluate
    for ( size_t i = 1 ; i < ncheck ; ++i ) {
      if ( res[i] > res[i-1] or res[i] > norm ) {
        stop = false;
        break;
      }
    }
    // - stop
    if ( stop ) return iiter;
  }
}

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
