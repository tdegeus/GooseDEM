/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/FrictionDEM

================================================================================================= */

#ifndef FRICTIONFEM_H
#define FRICTIONFEM_H

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

#define FRICTIONDEM_WORLD_VERSION 0
#define FRICTIONDEM_MAJOR_VERSION 0
#define FRICTIONDEM_MINOR_VERSION 1

#define FRICTIONDEM_VERSION_AT_LEAST(x,y,z) \
  (FRICTIONDEM_WORLD_VERSION>x || (FRICTIONDEM_WORLD_VERSION>=x && \
  (FRICTIONDEM_MAJOR_VERSION>y || (FRICTIONDEM_MAJOR_VERSION>=y && \
                              FRICTIONDEM_MINOR_VERSION>=z))))

#define FRICTIONDEM_VERSION(x,y,z) \
  (FRICTIONDEM_WORLD_VERSION==x && \
   FRICTIONDEM_MAJOR_VERSION==y && \
   FRICTIONDEM_MINOR_VERSION==z)

// ------------------------------------------ alias types ------------------------------------------

namespace FrictionDEM {

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
  typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic,              1, Eigen::ColMajor> ColS;

}

// -------------------------------------------------------------------------------------------------

#include "Constitutive.h"
#include "Geometry.h"
#include "Sim.h"

#include "Constitutive.cpp"
#include "Geometry.cpp"
#include "Sim.cpp"

// -------------------------------------------------------------------------------------------------

#endif
