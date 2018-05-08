#ifndef __GEMINO_RATE_MATRICES_H
#define __GEMINO_RATE_MATRICES_H



/*
 ******************************************************************************
 ******************************************************************************

 Create rate matrices for multiple molecular-evolution models:
 TN93 and its special cases (JC69, K80, F81, HKY85, and F84), plus GTR and UNREST

 ******************************************************************************
 ******************************************************************************
 */



#include <RcppArmadillo.h>

#include "gemino_types.h" // integer types

using namespace Rcpp;


arma::mat TN93_rate_matrix(const std::vector<double>& pi_tcag,
                           const double& alpha_1, const double& alpha_2,
                           const double& beta, const double& xi);

arma::mat JC69_rate_matrix(const double& lambda, const double& xi);

arma::mat K80_rate_matrix(const double& alpha, const double& beta,
                          const double& xi);

arma::mat F81_rate_matrix(const std::vector<double>& pi_tcag,
                          const double& xi);

arma::mat HKY85_rate_matrix(const std::vector<double>& pi_tcag,
                            const double& alpha, const double& beta,
                            const double& xi);

arma::mat F84_rate_matrix(const std::vector<double>& pi_tcag,
                          const double& beta, const double& kappa,
                          const double& xi);

arma::mat GTR_rate_matrix(const std::vector<double>& pi_tcag,
                          const std::vector<double>& abcdef,
                          const double& xi);

void UNREST_rate_matrix(arma::mat& Q, std::vector<double>& pi_tcag, const double& xi);



#endif
