
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
arma::mat TN93_rate_matrix(const std::vector<double>& pi_tcag,
                           const double& alpha_1, const double& alpha_2,
                           const double& beta, const double& xi) {

    arma::mat Q(4,4);
    Q.fill(beta);
    Q.submat(arma::span(0,1), arma::span(0,1)).fill(alpha_1);
    Q.submat(arma::span(2,3), arma::span(2,3)).fill(alpha_2);
    for (uint i = 0; i < 4; i++) Q.col(i) *= pi_tcag[i];

    // Filling in diagonals
    Q.diag().fill(0.0);  // reset to zero so summing by row works
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi * 0.25;
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
arma::mat JC69_rate_matrix(const double& lambda, const double& xi) {

    std::vector<double> pi_tcag = {1, 1, 1, 1};

    arma::mat Q = TN93_rate_matrix(pi_tcag, lambda, lambda, lambda, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the K80 substitution model.
//'
//' K80 is a special case of TN93.
//'
//' @noRd
//'
arma::mat K80_rate_matrix(const double& alpha, const double& beta,
                          const double& xi) {

    std::vector<double> pi_tcag = {1, 1, 1, 1};

    arma::mat Q = TN93_rate_matrix(pi_tcag, alpha, alpha, beta, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the F81 substitution model.
//'
//' F81 is a special case of TN93.
//'
//' @noRd
//'
arma::mat F81_rate_matrix(const std::vector<double>& pi_tcag, const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_tcag, 1, 1, 1, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the HKY85 substitution model.
//'
//' HKY85 is a special case of TN93.
//'
//' @noRd
//'
arma::mat HKY85_rate_matrix(const std::vector<double>& pi_tcag,
                            const double& alpha, const double& beta,
                            const double& xi) {

    arma::mat Q = TN93_rate_matrix(pi_tcag, alpha, alpha, beta, xi);

    return Q;
}


//' Q matrix for rates for a given nucleotide using the F84 substitution model.
//'
//' F84 is a special case of TN93.
//'
//' @noRd
//'
arma::mat F84_rate_matrix(const std::vector<double>& pi_tcag,
                          const double& beta, const double& kappa,
                          const double& xi) {

    double pi_y = pi_tcag[0] + pi_tcag[1];
    double pi_r = pi_tcag[2] + pi_tcag[3];

    double alpha_1 = 1 + kappa / pi_y;
    double alpha_2 = 1 + kappa / pi_r;

    arma::mat Q = TN93_rate_matrix(pi_tcag, alpha_1, alpha_2, beta, xi);

    return Q;
}




//' Q matrix for rates for a given nucleotide using the GTR substitution model.
//'
//' @noRd
//'
arma::mat GTR_rate_matrix(const std::vector<double>& pi_tcag,
                          const std::vector<double>& abcdef,
                          const double& xi) {

    arma::mat Q(4, 4, arma::fill::zeros);

    // Filling in non-diagonals
    uint k = 0;
    for (uint i = 0; i < 3; i++) {
        for (uint j = i+1; j < 4; j++) {
            Q(i,j) = abcdef[k];
            Q(j,i) = abcdef[k];
            k++;
        }
    }
    for (uint i = 0; i < 4; i++) Q.col(i) *= pi_tcag[i];

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums += xi * 0.25;
    rowsums *= -1;
    Q.diag() = rowsums;

    return Q;

}







//' Estimates equilibrium nucleotide frequencies from an input rate matrix.
//'
//' It does this by solving for πQ = 0 by finding the left eigenvector of Q that
//' corresponds to the eigenvalue closest to zero.
//' This is only needed for the UNREST model.
//'
//' @inheritParams Q UNREST_rate_matrix
//' @inheritParams pi_tcag UNREST_rate_matrix
//'
//' @noRd
//'
inline void est_pi_tcag(const arma::mat& Q, std::vector<double>& pi_tcag) {

    arma::cx_vec eigvals;
    arma::cx_mat eigvecs;

    arma::eig_gen(eigvals, eigvecs, Q.t());

    arma::vec vals = arma::abs(arma::real(eigvals));
    arma::mat vecs = arma::real(eigvecs);

    uint i = arma::as_scalar(arma::find(vals == arma::min(vals), 1));

    arma::vec left_vec = vecs.col(i);
    double sumlv = arma::accu(left_vec);

    pi_tcag.resize(4);
    for (uint i = 0; i < 4; i++) pi_tcag[i] = left_vec(i) / sumlv;

    return;
}





//' Q matrix for rates for a given nucleotide using the UNREST substitution model.
//'
//' This function also fills in a vector of equilibrium frequencies for each nucleotide.
//' This calculation has to be done for this model only because it uses separate
//' values for each non-diagonal cell and doesn't use equilibrium frequencies for
//' creating the matrix.
//'
//'
//' @param Q Matrix of substitution rates for "T", "C", "A", and "G", respectively.
//'     Do not include indel rates here! Diagonal values are ignored.
//' @param pi_tcag Empty vector of equilibrium frequencies for for "T", "C", "A", and "G",
//'     respectively. This vector will be filled in by this function.
//' @param xi Overall rate of indels.
//'
//' @noRd
//'
void UNREST_rate_matrix(arma::mat& Q, std::vector<double>& pi_tcag, const double& xi) {

    if (Q.n_rows != 4 || Q.n_cols != 4) stop("Q matrix should be 4 x 4");

    // reset to zero so summing by row works
    Q.diag().fill(0.0);
    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums *= -1;
    Q.diag() = rowsums;

    // Estimate pi_tcag before incorporating indels:
    est_pi_tcag(Q, pi_tcag);

    /*
     Now include indel rates
     (I'm subtracting bc diagonals are negative)
     (Also, I'm doing this here so it doesn't affect pi_tcag estimates)
     */
    Q.diag() -= xi * 0.25;

    return;
}


