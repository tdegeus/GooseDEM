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

inline void velocityVerlet(Geometry &g, double dt, const ColS &iip, const ColD &vp, const ColD &ap)
{
  // local variable and history
  ColD V;
  ColD V_n  = g.dofs_v();
  ColD A;
  ColD A_n  = g.dofs_a();
  ColD M    = g.dofs_m();
  ColD Minv = M.cwiseInverse();

  // (1) new positions: x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  g.set_x( g.x() + dt * g.v() + 0.5 * dt * g.a() );

  // (2a) estimate new velocities
  // - update velocities (DOFs)
  V = g.dofs_v() + dt * g.dofs_a();
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = vp(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - update velocities (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = vp(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );

  // (2b) new velocities
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - update velocities (DOFs)
  V = V_n + .5 * dt * ( A_n + A );
  // - set prescribed velocities
  for ( int i = 0 ; i < iip.size() ; ++i ) V(iip(i)) = vp(i);
  // - update velocities (particles)
  g.set_v( g.asParticle(V) );

  // (3) new accelerations
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - g.dofs_f() );
  // - set prescribed accelerations
  for ( int i = 0 ; i < iip.size() ; ++i ) A(iip(i)) = ap(i);
  // - update accelerations (particles)
  g.set_a( g.asParticle(A) );
}

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
