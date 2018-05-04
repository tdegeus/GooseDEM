/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_EXT_FRICTION_H
#define GOOSEDEM_EXT_FRICTION_H

// -------------------------------------------------------------------------------------------------

#include "../../GooseDEM.h"

// -------------------------------------- version information --------------------------------------

#define GOOSEDEM_EXT_FRICTION_WORLD_VERSION 0
#define GOOSEDEM_EXT_FRICTION_MAJOR_VERSION 0
#define GOOSEDEM_EXT_FRICTION_MINOR_VERSION 1

#define GOOSEDEM_EXT_FRICTION_VERSION_AT_LEAST(x,y,z) \
  (GOOSEDEM_EXT_FRICTION_WORLD_VERSION>x || (GOOSEDEM_EXT_FRICTION_WORLD_VERSION>=x && \
  (GOOSEDEM_EXT_FRICTION_MAJOR_VERSION>y || (GOOSEDEM_EXT_FRICTION_MAJOR_VERSION>=y && \
                                GOOSEDEM_EXT_FRICTION_MINOR_VERSION>=z))))

#define GOOSEDEM_EXT_FRICTION_VERSION(x,y,z) \
  (GOOSEDEM_EXT_FRICTION_WORLD_VERSION==x && \
   GOOSEDEM_EXT_FRICTION_MAJOR_VERSION==y && \
   GOOSEDEM_EXT_FRICTION_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#include "PotentialAdhesion.h"
#include "PotentialAdhesion.cpp"

#include "Geometry.h"
#include "Geometry.cpp"

// -------------------------------------------------------------------------------------------------

#endif
