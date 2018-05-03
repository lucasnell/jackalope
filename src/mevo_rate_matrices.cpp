
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
#include "mevo_rate_matrices.h"

using namespace Rcpp;



//' Q matrix for rates for a given nucleotide using the TN93 substitution model.
//'
//' @noRd
//'
inline arma::mat TN93_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& alpha_1, const double& alpha_2, const double& beta,
        const double& xi) {

    arma::vec pis = {pi_t, pi_c, pi_a, pi_g};

    arma::mat Q(4,4);
    Q.fill(beta);
    Q.submat(arma::span(0,1), arma::span(0,1)).fill(alpha_1);
    Q.submat(arma::span(2,3), arma::span(2,3)).fill(alpha_2);
    for (uint i = 0; i < 4; i++) Q.col(i) *= pis(i);

    // Filling in diagonals
    Q.diag().fill(0.0);  // reset to zero so summing by row works
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;
}



//' Q matrix for rates for a given nucleotide using the JC69 substitution model.
//'
//' JC69 is a special case of TN93.
//'
//' @noRd
//'
inline arma::mat JC69_rate_matrix(
        const double& lambda, const double& xi) {

    arma::mat Q = TN93_rate_matrix(1, 1, 1, 1, lambda, lambda, lambda, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the K80 substitution model.
//'
//' K80 is a special case of TN93.
//'
//' @noRd
//'
inline arma::mat K80_rate_matrix(
        const double& alpha, const double& beta,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(1, 1, 1, 1, alpha, alpha, beta, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the F81 substitution model.
//'
//' F81 is a special case of TN93.
//'
//' @noRd
//'
inline arma::mat F81_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, 1, 1, 1, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the HKY85 substitution model.
//'
//' HKY85 is a special case of TN93.
//'
//' @noRd
//'
inline arma::mat HKY85_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& alpha, const double& beta,
        const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha, alpha, beta, xi);


    return Q;
}


//' Q matrix for rates for a given nucleotide using the F84 substitution model.
//'
//' F84 is a special case of TN93.
//'
//' @noRd
//'
inline arma::mat F84_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& beta, const double& kappa,
        const double& xi) {

    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    double alpha_1 = 1 + kappa / pi_y;
    double alpha_2 = 1 + kappa / pi_r;

    arma::mat Q = TN93_rate_matrix(pi_a, pi_c, pi_g, pi_t, alpha_1, alpha_2, beta, xi);

    return Q;
}




//' Q matrix for rates for a given nucleotide using the GTR substitution model.
//'
//' @noRd
//'
inline arma::mat GTR_rate_matrix(
        const double& pi_t, const double& pi_c,
        const double& pi_a, const double& pi_g,
        const double& a, const double& b, const double& c,
        const double& d, const double& e, const double& f,
        const double& xi) {

    arma::vec pis = {pi_t, pi_c, pi_a, pi_g};

    arma::mat Q(4, 4, arma::fill::zeros);

    // Filling in non-diagonals
    arma::vec letters = {a, b, c, d, e, f};
    uint k = 0;
    for (uint i = 0; i < 3; i++) {
        for (uint j = i+1; j < 4; j++) {
            Q(i,j) = letters(k);
            Q(j,i) = letters(k);
            k++;
        }
    }
    for (uint i = 0; i < 4; i++) Q.col(i) *= pis(i);

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;

}


//' Q matrix for rates for a given nucleotide using the UNREST substitution model.
//'
//' @param Q Matrix of rates for "T", "C", "A", and "G", respectively.
//'     Diagonal values are ignored.
//' @param xi Overall rate of indels.
//'
//' @noRd
//'
inline arma::mat UNREST_rate_matrix(
        arma::mat Q, const double& xi) {

    // reset to zero so summing by row works
    Q.diag().fill(0.0);
    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;

}


