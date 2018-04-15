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

  MatS m_particles; // 'list' particle-id pairs   [n, ndim]
  ColD m_k;         // stiffness                  [n]
  ColD m_D0;        // relaxed length             [n]

public:

  // constructor
  Spring(){};
  Spring(const MatS &particles, const ColD &k, const ColD &D0);

  // compute the force on each particle
  MatD force(const MatD &x) const;

  // compute the coordination
  ColS coordination(const MatD &X) const;

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
