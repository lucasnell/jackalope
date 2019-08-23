#ifndef __JACKALOPE_TYPES_H
#define __JACKALOPE_TYPES_H

#define __JACKALOPE_DEBUG // comment out when done debugging



/*
 ********************************************************

 Basic integer types used throughout

 ********************************************************
 */

#ifndef __JACKALOPE_DEBUG
#define ARMA_NO_DEBUG
#endif

#include <RcppArmadillo.h>
#include <cstdint>
#include <vector>
#include <string>
#include <deque>
#include "pcg/pcg_extras.hpp"  // pcg 128-bit integer type




using namespace Rcpp;

typedef uint_fast8_t uint8;
typedef uint_fast32_t uint32;
typedef int_fast32_t sint32;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;
typedef pcg_extras::pcg128_t uint128;



namespace jlp {
const std::string bases = "TCAG";
}



#endif
