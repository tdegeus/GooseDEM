/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_GEOMETRY_CPP
#define GOOSEDEM_GEOMETRY_CPP

// -------------------------------------------------------------------------------------------------

#include "Geometry.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

inline Geometry::Geometry(const ColD &m, const MatD &x, const MatS &dofs) :
  m_x(x), m_m(m), m_dofs(dofs)
{
  // extract dimensions
  m_N    = static_cast<size_t>(m_dofs.rows());
  m_ndim = static_cast<size_t>(m_dofs.cols());
  m_ndof = static_cast<size_t>(m_dofs.maxCoeff() + 1);

  // check input
  assert( m_x.rows() == m_N    );
  assert( m_x.cols() == m_ndim );
  assert( m_m.size() == m_N    );

  // allocate particle vectors
  m_v = MatD(m_x.rows(), m_x.cols());
  m_a = MatD(m_x.rows(), m_x.cols());

  // zero-initialize particle vectors
  m_v.setZero();
  m_a.setZero();
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set(const Spring &mat)
{
  m_spring = mat;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set(const Dashpot &mat)
{
  m_dashpot = mat;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::x() const
{
  return m_x;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::v() const
{
  return m_v;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::a() const
{
  return m_a;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::m() const
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::f() const
{
  // allocate force vector per particle
  MatD f(m_N, m_ndim);

  // zero-initialize
  f.setZero();

  // evaluate constitutive models
  f += m_spring .force(m_x);
  f += m_dashpot.force(m_x);

  return f;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_v() const
{
  return asDofs(m_v);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_a() const
{
  return asDofs(m_a);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_m() const
{
  return asDofs(m_m);
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::dofs_f() const
{
  return assembleDofs(f());
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_x(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_x = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_v(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_v = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_a(const MatD &pvector)
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // store
  m_a = pvector;
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_v(const ColD &dofval)
{
  set_v(asParticle(dofval));
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::set_a(const ColD &dofval)
{
  set_a(asParticle(dofval));
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::asDofs(const ColD &pscalar) const
{
  // check input
  assert( static_cast<size_t>(pscalar.rows()) == m_N );

  // list of DOF values
  // - allocate
  ColD dofval(m_ndof);
  // - zero-initialize
  dofval.setZero();

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = pscalar(n);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::asDofs(const MatD &pvector) const
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // list of DOF values
  // - allocate
  ColD dofval(m_ndof);
  // - zero-initialize
  dofval.setZero();

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = pvector(n,i);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline ColD Geometry::assembleDofs(const MatD &pvector) const
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // list of DOF values
  // - allocate
  ColD dofval(m_ndof);
  // - zero-initialize
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; list of DOF values
    // - allocate
    ColD t_dofval(m_ndof);
    // - zero-initialize
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_N ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        t_dofval(m_dofs(n,i)) += pvector(n,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline MatD Geometry::asParticle(const ColD &dofval) const
{
  // check input
  assert( dofval.size() == m_ndof );

  // list of particle vectors
  // - allocate
  MatD pvector(m_N, m_ndim);
  // - zero-initialize
  pvector.setZero();

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      pvector(n,i) = dofval(m_dofs(n,i));

  return pvector;
}

// -------------------------------------------------------------------------------------------------

inline void dump(const std::string &fname, const ColD &matrix)
{
  const static Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", "\n");

  std::ofstream file(fname.c_str());

  file << matrix.format(fmt);
}

// -------------------------------------------------------------------------------------------------

inline void dump(const std::string &fname, const MatD &matrix)
{
  const static Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", "\n");

  std::ofstream file(fname.c_str());

  file << matrix.format(fmt);
}

// -------------------------------------------------------------------------------------------------

}

// -------------------------------------------------------------------------------------------------

#endif
