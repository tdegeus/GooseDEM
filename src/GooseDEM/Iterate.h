/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_ITERATE_H
#define GOOSEDEM_ITERATE_H

// -------------------------------------------------------------------------------------------------

#include "GooseDEM.h"

// =========================================== GooseDEM ============================================

namespace GooseDEM {

// -------------------------------------------------------------------------------------------------

class StopList
{
private:

  // list with residuals
  std::vector<double> m_res;

public:

  // constructors
  StopList(size_t n=1);

  // reset all residuals to infinity (and change the number of residuals to check)
  void reset();
  void reset(size_t n);

  // update list of residuals, return true if all residuals are below the tolerance
  bool stop(double res, double tol);

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif
