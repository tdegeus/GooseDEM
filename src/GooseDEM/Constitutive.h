/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_CONSTITUTIVE_H
#define GOOSEDEM_CONSTITUTIVE_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Spring
{
private:

  MatS m_particles; // 'list' particle-id pairs   [nspring, ndim]
  ColD m_k;         // stiffness                  [nspring]
  ColD m_D0;        // relaxed length             [nspring]

public:

  // constructor
  Spring(){};
  Spring(const MatS &particles, const ColD &k, const ColD &D0);

  // compute the spring force on each particle
  MatD force(const MatD &x);

  // compute the coordination number of each spring
  ColS coordination(const MatD &X);

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
