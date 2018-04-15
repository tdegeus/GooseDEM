/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseDEM

================================================================================================= */

#ifndef GOOSEDEM_H
#define GOOSEDEM_H

// --------------------------------------- include libraries ---------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

// ---------------------------------------- dummy operation ----------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------- version information --------------------------------------

#define GOOSEDEM_WORLD_VERSION 0
#define GOOSEDEM_MAJOR_VERSION 0
#define GOOSEDEM_MINOR_VERSION 1

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

#include "Spring.h"
#include "Geometry.h"
#include "Sim.h"

#include "Spring.cpp"
#include "Geometry.cpp"
#include "Sim.cpp"

// -------------------------------------------------------------------------------------------------

#endif
