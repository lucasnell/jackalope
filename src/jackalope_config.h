#ifndef __JACKALOPE_CONFIG_H
#define __JACKALOPE_CONFIG_H

/*
 Include this file in all C++ files to control debugging and verbose output for
 diagnostics.
 It must go before any other includes (especially RcppArmadillo) so it can control
 Armadillo's debugging.
 */


// Comment this out when done debugging:
#define __JACKALOPE_DEBUG
// Comment this out when done with diagnostics:
#define __JACKALOPE_DIAGNOSTICS


#ifndef __JACKALOPE_DEBUG
#define ARMA_NO_DEBUG
#endif



#endif

