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
  Geometry &geometry, double dt, const ColS &iip=ColS(), const ColD &v_p=ColD());

// iterate until all particles have come to a rest
inline size_t quasiStaticVelocityVerlet(
  Geometry &g, double dt, const ColS &iip, const ColD &v_p, double norm, size_t ncheck=20)

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
