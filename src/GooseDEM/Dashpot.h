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

  MatS m_particles; // particle pairs   [n, ndim]
  ColD m_eta;       // damping constant [n]

public:

  // constructor
  Dashpot(){};
  Dashpot(const MatS &particles, const ColD &eta);

  // compute the force on each particle (the output could contain many zero rows)
  MatD force(const MatD &v) const;

  // compute the coordination of each particle
  ColS coordination(const MatD &X) const;

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
