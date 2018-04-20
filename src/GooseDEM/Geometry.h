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
private:

  // particles
  MatD m_x;     // position     [N, ndim]
  MatD m_v;     // velocity     [N, ndim]
  MatD m_a;     // acceleration [N, ndim]
  ColD m_m;     // mass         [N, ndim]
  MatS m_dofs;  // DOF-number   [N, ndim]

  // dimensions
  size_t m_N;    // number of particles
  size_t m_ndim; // number of spatial dimensions
  size_t m_ndof; // number of DOFs

  // constitutive models
  Spring  m_spring;
  Dashpot m_dashpot;

public:

  // constructor
  Geometry(const ColD &m, const MatD &x, const MatS &dofs);

  // append constitutive models
  void set(const Spring  &mat);
  void set(const Dashpot &mat);

  // return particle vectors [N, ndim] or values [N]
  MatD x() const;
  MatD v() const;
  MatD a() const;
  MatD f() const;
  ColD m() const;

  // return DOF values [ndof]
  ColD dofs_v() const;
  ColD dofs_a() const;
  ColD dofs_f() const;
  ColD dofs_m() const;

  // overwrite particle vectors [N, ndim]
  void set_x(const MatD &pvector);
  void set_v(const MatD &pvector);
  void set_a(const MatD &pvector);

  // overwrite particle vectors, reconstructed from DOF values
  void set_v(const ColD &dofval); // == set_v(asParticle(V))
  void set_a(const ColD &dofval); // == set_a(asParticle(A))

  // convert to DOF values (overwrite entries that occur more that once): [N, ndim] -> [ndof]
  ColD asDofs(const ColD &pscalar) const;
  ColD asDofs(const MatD &pvector) const;

  // assemble vectors (adds entries that occur more that once): [N, ndim] -> [ndof]
  ColD assembleDofs(const MatD &pvector) const;

  // reconstruct particle vectors: [ndof] -> [N, ndim]
  MatD asParticle(const ColD &dofval) const;

};

// -------------------------------------------------------------------------------------------------

inline void dump(const std::string &fname, const ColD &matrix);
inline void dump(const std::string &fname, const MatD &matrix);

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
