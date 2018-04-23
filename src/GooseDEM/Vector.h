/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_VECTOR_H
#define GOOSEDEM_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// =========================================== GooseDEM ============================================

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Vector
{
private:

  // particles
  MatS m_dofs;   // DOF-number [N, ndim]

  // dimensions
  size_t m_N;    // number of particles
  size_t m_ndim; // number of spatial dimensions
  size_t m_ndof; // number of DOFs

public:

  // constructor
  Vector(){};
  Vector(const MatS &dofs);

  // convert to DOF values (overwrite entries that occur more that once): [N, ndim] -> [ndof]
  ColD asDofs(const ColD &pscalar) const;
  ColD asDofs(const MatD &pvector) const;

  // reconstruct particle vectors: [ndof] -> [N, ndim]
  MatD asParticle(const ColD &dofval) const;

  // assemble vectors (adds entries that occur more that once): [N, ndim] -> [ndof]
  ColD assembleDofs(const MatD &pvector) const;

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif
