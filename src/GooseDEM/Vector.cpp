/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_VECTOR_CPP
#define GOOSEDEM_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "Vector.h"

// =========================================== GooseDEM ============================================

namespace GooseDEM {

// ------------------------------------------ constructor ------------------------------------------

inline Vector::Vector(const MatS &dofs) : m_dofs(dofs)
{
  // extract dimensions
  m_N    = static_cast<size_t>(m_dofs.rows());
  m_ndim = static_cast<size_t>(m_dofs.cols());
  m_ndof = static_cast<size_t>(m_dofs.maxCoeff() + 1);
}

// ----------------------------------- particle scalar -> dofval -----------------------------------

inline ColD Vector::asDofs(const ColD &pscalar) const
{
  // check input
  assert( static_cast<size_t>(pscalar.rows()) == m_N );

  // list of DOF values
  ColD dofval = ColD::Zero(m_ndof);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = pscalar(n);

  return dofval;
}

// ----------------------------------- particle vector -> dofval -----------------------------------

inline ColD Vector::asDofs(const MatD &pvector) const
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // list of DOF values
  ColD dofval = ColD::Zero(m_ndof);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = pvector(n,i);

  return dofval;
}

// ----------------------------------- dofval -> particle vector -----------------------------------

inline MatD Vector::asParticle(const ColD &dofval) const
{
  // check input
  assert( dofval.size() == m_ndof );

  // list of particle vectors
  MatD pvector = MatD::Zero(m_N, m_ndim);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_N ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      pvector(n,i) = dofval(m_dofs(n,i));

  return pvector;
}

// ----------------------------------- particle vector -> dofval -----------------------------------

inline ColD Vector::assembleDofs(const MatD &pvector) const
{
  // check input
  assert( static_cast<size_t>(pvector.rows()) == m_N    );
  assert( static_cast<size_t>(pvector.cols()) == m_ndim );

  // list of DOF values
  ColD dofval = ColD::Zero(m_ndof);

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // list of DOF values
    ColD t_dofval = ColD::Zero(m_ndof);

    // assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_N ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        t_dofval(m_dofs(n,i)) += pvector(n,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif
