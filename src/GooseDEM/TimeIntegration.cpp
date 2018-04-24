/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_TIMEINTEGRATION_CPP
#define GOOSEDEM_TIMEINTEGRATION_CPP

// -------------------------------------------------------------------------------------------------

#include "TimeIntegration.h"

// =================================================================================================

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline void Verlet(Geometry &g, double dt)
{
  // history

  ColD V;
  ColD A;
  ColD V_n = g.dofs_v();
  ColD A_n = g.dofs_a();

  // new position

  g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

  // new velocity

  A = g.solve();

  g.set_a( A );

  // new velocity

  g.set_v( V_n + .5 * dt * ( A_n + A ) );

  // finalize

  g.timestep(dt);
}

// -------------------------------------------------------------------------------------------------

inline void velocityVerlet(Geometry &g, double dt)
{
  // history

  ColD V;
  ColD A;
  ColD V_n = g.dofs_v();
  ColD A_n = g.dofs_a();

  // new position

  g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

  // estimate new velocity

  g.set_v( V_n + dt * A_n );

  A = g.solve();

  g.set_v( V_n + .5 * dt * ( A_n + A ) );

  // new velocity

  A = g.solve();

  g.set_v( V_n + .5 * dt * ( A_n + A ) );

  // new acceleration

  A = g.solve();

  g.set_a(A);

  // finalize

  g.timestep(dt);
}

// -------------------------------------------------------------------------------------------------

inline size_t quasiStaticVelocityVerlet(Geometry &g, double dt, double tol)
{
  // reset residuals
  g.reset();

  // zero-initialize iteration counter
  size_t iiter = 0;

  // loop until convergence
  while ( true )
  {
    // - update iteration counter
    iiter++;
    // - time-step
    velocityVerlet(g,dt);
    // - check for convergence
    if ( g.stop(tol) ) return iiter;
  }
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif
