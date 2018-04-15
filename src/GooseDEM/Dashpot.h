/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_DASHPOT_H
#define GOOSEDEM_DASHPOT_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Dashpot
{
private:

  MatS m_particles; // 'list' particle-id pairs   [n, ndim]
  ColD m_eta;       // damping constant           [n]

public:

  // constructor
  Dashpot(){};
  Dashpot(const MatS &particles, const ColD &eta);

  // compute the force on each particle
  MatD force(const MatD &v) const;

  // compute the coordination
  ColS coordination(const MatD &X) const;

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
