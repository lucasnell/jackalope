# ifndef __JACKAL_TYPES_H
# define __JACKAL_TYPES_H


/*
 ********************************************************

 Basic integer types used throughout

 ********************************************************
 */



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




#endif
