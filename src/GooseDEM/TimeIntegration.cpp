/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_TIMEINTEGRATION_CPP
#define GOOSEDEM_TIMEINTEGRATION_CPP

// -------------------------------------------------------------------------------------------------

#include "TimeIntegration.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline void velocityVerlet(Geometry &g, double dt)
{
  // local variables and history
  ColD V;
  ColD A;
  ColD V_n = g.dofs_v();
  ColD A_n = g.dofs_a();

  // (1) new positions
  g.set_x( g.x() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() );

  // (2a) estimate new velocity
  // - new velocity
  V = V_n + dt * A_n;
  // - update particles
  g.set_v(V);
  // - solve for accelerations
  A = g.solve();
  // - new velocity
  V = V_n + .5 * dt * ( A_n + A );
  // - update particles
  g.set_v(V);

  // (2b) new velocity
  // - solve for accelerations
  A = g.solve();
  // - new velocity
  V = V_n + .5 * dt * ( A_n + A );
  // - update particles
  g.set_v(V);

  // (3) new accelerations
  // - solve for accelerations
  A = g.solve();
  // - update particles
  g.set_a(A);

  // process time-step
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

}

// -------------------------------------------------------------------------------------------------

#endif
