/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_H
#define GOOSEDEM_H

// --------------------------------------- include libraries ---------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iso646.h>
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

// ---------------------------------------- dummy operation ----------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------- version information --------------------------------------

#define GOOSEDEM_WORLD_VERSION 0
#define GOOSEDEM_MAJOR_VERSION 0
#define GOOSEDEM_MINOR_VERSION 3

#define GOOSEDEM_VERSION_AT_LEAST(x,y,z) \
  (GOOSEDEM_WORLD_VERSION>x || (GOOSEDEM_WORLD_VERSION>=x && \
  (GOOSEDEM_MAJOR_VERSION>y || (GOOSEDEM_MAJOR_VERSION>=y && \
                                GOOSEDEM_MINOR_VERSION>=z))))

#define GOOSEDEM_VERSION(x,y,z) \
  (GOOSEDEM_WORLD_VERSION==x && \
   GOOSEDEM_MAJOR_VERSION==y && \
   GOOSEDEM_MINOR_VERSION==z)

// ------------------------------------------ alias types ------------------------------------------

namespace GooseDEM {

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
  typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic,              1, Eigen::ColMajor> ColS;

}

// -------------------------------------------------------------------------------------------------

#include "Write.h"
#include "Spring.h"
#include "Dashpot.h"
#include "Iterate.h"
#include "Vector.h"
#include "Geometry.h"
#include "TimeIntegration.h"

#include "Write.cpp"
#include "Spring.cpp"
#include "Dashpot.cpp"
#include "Iterate.cpp"
#include "Vector.cpp"
#include "TimeIntegration.cpp"

// -------------------------------------------------------------------------------------------------

#endif
