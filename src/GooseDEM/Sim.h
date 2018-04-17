/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_SIM_H
#define GOOSEDEM_SIM_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

// evaluate one time step
inline void velocityVerlet(
  Geometry &geometry, double dt,
  const ColS &ivp=ColS(), const ColD &vp=ColD(), const ColS &ifp=ColS(), const ColD &fp=ColD()
);

// iterate until all particles have come to a rest
inline size_t quasiStaticVelocityVerlet(
  Geometry &geometry, double dt, double norm, size_t ncheck=20,
  const ColS &ivp=ColS(), const ColD &vp=ColD(), const ColS &ifp=ColS(), const ColD &fp=ColD()
);

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
