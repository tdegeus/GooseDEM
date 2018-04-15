/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_SIM_H
#define GOOSEDEM_SIM_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class Geometry
{
private:

  // particles
  MatD m_x;     // position     [N, ndim]
  MatD m_v;     // velocity     [N, ndim]
  MatD m_a;     // acceleration [N, ndim]
  ColD m_m;     // mass         [N, ndim]
  MatS m_dofs;  // DOF-number   [N, ndim]

  // constitutive models
  Spring m_spring;

public:

  // constructor
  Geometry(const ColD &m, const MatD &x, const MatS &dofs);

  // append constitutive models
  void setSpring(const Spring &mat);

  // return particle values
  MatD x(); // position vector     [N, ndim]
  MatD v(); // velocity vector     [N, ndim]
  MatD a(); // acceleration vector [N, ndim]
  MatD f(); // force vector        [N, ndim]
  ColD m(); // mass                [N, ndim]

  // return DOF values
  ColD dofs_x(); // position     [ndof]
  ColD dofs_v(); // velocity     [ndof]
  ColD dofs_a(); // acceleration [ndof]
  ColD dofs_f(); // force        [ndof]
  ColD dofs_m(); // mass         [ndof]

  // overwrite particle values
  void set_x(const MatD &pvec); // position vector     [N, ndim]
  void set_v(const MatD &pvec); // velocity vector     [N, ndim]
  void set_a(const MatD &pvec); // acceleration vector [N, ndim]

  // reconstruct particle values from DOF values
  void set_dofs_x(const MatD &pvec); // position     [ndof]
  void set_dofs_v(const MatD &pvec); // velocity     [ndof]
  void set_dofs_a(const MatD &pvec); // acceleration [ndof]

};

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
