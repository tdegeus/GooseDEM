/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_SPRING_H
#define GOOSEDEM_SPRING_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Spring
{
private:

  MatS m_particles; // particle pairs   [n, ndim]
  ColD m_k;         // stiffness        [n]
  ColD m_D0;        // relaxed length   [n]

public:

  // constructor
  Spring(){};
  Spring(const MatS &particles, const ColD &k, const ColD &D0);

  // compute the force on each particle (the output could contain many zero rows)
  MatD force(const MatD &x) const;

  // compute the coordination of each particle
  ColS coordination(const MatD &X) const;

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
