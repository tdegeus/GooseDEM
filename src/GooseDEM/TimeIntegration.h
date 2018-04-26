/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_TIMEINTEGRATION_H
#define GOOSEDEM_TIMEINTEGRATION_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// =================================================================================================

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

// evaluate one time step
inline void Verlet(GeometryFriction &geometry, double dt);

// evaluate one time step
inline void velocityVerlet(GeometryFriction &geometry, double dt);

// iterate until all particles have come to a rest
inline size_t quasiStaticVelocityVerlet(GeometryFriction &geometry, double dt, double tol);

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif
