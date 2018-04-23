/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_WRITE_CPP
#define GOOSEDEM_WRITE_CPP

// -------------------------------------------------------------------------------------------------

#include "Write.h"

// -------------------------------------------------------------------------------------------------

namespace GooseDEM {

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
