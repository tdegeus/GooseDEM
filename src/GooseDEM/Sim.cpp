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

inline void velocityVerlet(Geometry &g, double dt,
  const ColS &ivp, const ColD &vp, const ColS &ifp, const ColD &fp)
{
  // local variables and history
  ColD F, V, A;
  ColD V_n  = g.dofs_v();
  ColD A_n  = g.dofs_a();
  ColD M    = g.dofs_m();
  ColD Minv = M.cwiseInverse();

  // (1) new positions: x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

  // (2a) estimate new velocity
  // - update velocity (DOFs)
  V = V_n + dt * A_n;
  // - set prescribed velocity
  for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
  // - update velocity (particles)
  g.set_v( g.asParticle(V) );
  // - get force
  F = g.dofs_f();
  // - set prescribed force
  for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( F );
  // - update velocity (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocity
  for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
  // - update velocity (particles)
  g.set_v( g.asParticle(V) );

  // (2b) new velocity
  // - get force
  F = g.dofs_f();
  // - set prescribed force
  for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( F );
  // - update velocity (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocity
  for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
  // - update velocity (particles)
  g.set_v( g.asParticle(V) );

  // (3) new accelerations
  // - get force
  F = g.dofs_f();
  // - set prescribed force
  for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( F );
  // - enforce prescribed velocity
  for ( auto i = 0 ; i < ivp.size() ; ++i ) A(ivp(i)) = 0.0;
  // - update accelerations (particles)
  g.set_a( g.asParticle(A) );
}

// -------------------------------------------------------------------------------------------------

inline size_t quasiStaticVelocityVerlet(
  Geometry &g, double dt, double norm, size_t ncheck,
  const ColS &ivp, const ColD &vp, const ColS &ifp, const ColD &fp)
{
  // local variables
  ColD   F, V, V_n, A, A_n;
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

    // (2a) estimate new velocity
    // - update velocity (DOFs)
    V = V_n + dt * A_n;
    // - set prescribed velocity
    for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
    // - update velocity (particles)
    g.set_v( g.asParticle(V) );
      // - get force
    F = g.dofs_f();
    // - set prescribed force
    for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( F );
    // - update velocity (DOFs)
    V = V_n + .5 * dt * ( A_n + A );
    // - set prescribed velocity
    for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
    // - update velocity (particles)
    g.set_v( g.asParticle(V) );

    // (2b) new velocity
      // - get force
    F = g.dofs_f();
    // - set prescribed force
    for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( F );
    // - update velocity (DOFs)
    V = V_n + .5 * dt * ( A_n + A );
    // - set prescribed velocity
    for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = vp(i);
    // - update velocity (particles)
    g.set_v( g.asParticle(V) );

    // (3) new accelerations
      // - get force
    F = g.dofs_f();
    // - set prescribed force
    for ( auto i = 0 ; i < ifp.size() ; ++i ) F(ifp(i)) += fp(i);
    // - solve for accelerations (DOFs)
    A = Minv.cwiseProduct( F );
    // - enforce prescribed velocity
    for ( auto i = 0 ; i < ivp.size() ; ++i ) A(ivp(i)) = 0.0;
    // - update accelerations (particles)
    g.set_a( g.asParticle(A) );

    // check for convergence
    // ---------------------

    // compute total kinetic energy
    // - ignore prescribed DOFs
    for ( auto i = 0 ; i < ivp.size() ; ++i ) V(ivp(i)) = 0.0;
    // - zero-initialize
    new_res = 0.0;
    // - compute
    for ( auto i = 0 ; i < V.size() ; ++i )
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
