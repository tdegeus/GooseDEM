/* =================================================================================================

Contributors:

  Tianxia Ma
  Wencheng Ji

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_EXT_FRICTION_POTENTIALADHESION_H
#define GOOSEDEM_EXT_FRICTION_POTENTIALADHESION_H

// -------------------------------------------------------------------------------------------------

#include "Friction.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {
namespace Ext {
namespace Friction {

// -------------------------------------------------------------------------------------------------

class PotentialAdhesion
{
private:

  MatS m_particles; // particle pairs   [n, ndim]
  ColD m_k;         // stiffness        [n]
  ColD m_D0;        // relaxed length   [n]

public:

  // constructor
  PotentialAdhesion(){};
  PotentialAdhesion(const MatS &particles, const ColD &k, const ColD &D0);

  // compute the force on each particle (the output could contain many zero rows)
  MatD force(const MatD &x) const;

  // compute the coordination of each particle
  ColS coordination(const MatD &X) const;

};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
