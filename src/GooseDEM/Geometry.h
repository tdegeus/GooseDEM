/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_GEOMETRY_H
#define GOOSEDEM_GEOMETRY_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Geometry
{
public:

  // solve for DOF-accelerations [ndof]
  virtual ColD solve() const { return ColD(); };

  // reset residuals, check for convergence
  virtual void reset()          { return;       }
  virtual bool stop(double tol) { return false; };

  // process time-step
  virtual void timestep(double dt) { return; };

  // return particle vectors [N, ndim]
  virtual MatD x() const { return MatD(); };
  virtual MatD v() const { return MatD(); };
  virtual MatD a() const { return MatD(); };

  // return DOF values [ndof]
  virtual ColD dofs_v() const { return ColD(); };
  virtual ColD dofs_a() const { return ColD(); };

  // overwrite particle vectors [N, ndim]
  virtual void set_x(const MatD &pvector) { return; };

  // overwrite particle vectors, reconstructed from DOF values [ndof]
  virtual void set_v(const ColD &dofval) { return; };
  virtual void set_a(const ColD &dofval) { return; };

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
