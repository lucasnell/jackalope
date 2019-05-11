
/*
 ******************************************************************************
 ******************************************************************************

 Create rate matrices for multiple molecular-evolution models:
 TN93 and its special cases (JC69, K80, F81, HKY85, and F84), plus GTR and UNREST

 ******************************************************************************
 ******************************************************************************
 */



#include <RcppArmadillo.h>

#include "jackalope_types.h" // integer types
#include "util.h" // str_stop

using namespace Rcpp;



/*
 Check that vectors meet criteria. Used for `pi_tcag` and `abcdef` below.
 */
inline void vec_check(const std::vector<double>& in_vec,
                      const std::string& vec_name,
                      const bool& zero_check,
                      const uint64& needed_size) {
    if (in_vec.size() != needed_size) {
        str_stop({"\nFor substitution models, the vector `", vec_name,
                 "` should always be of length ", std::to_string(needed_size), "."});
    }
    bool all_zero = true;
    for (const double& d : in_vec) {
        if (d < 0) {
            str_stop({"\nFor substitution models, all values in vector `", vec_name,
                     "` should be >= 0."});
        }
        if (d > 0) all_zero = false;
    }
    if (zero_check && all_zero) {
        str_stop({"\nFor substitution models, at least one value in vector `", vec_name,
                 "` should be > 0."});
    }
    return;
}



//' Construct necessary information for substitution models.
//'
//' For a more detailed explanation, see `vignette("sub-models")`.
//'
//'
//' @name sub_models
//'
//' @seealso \code{\link{create_variants}}
//'
//' @return A `sub_model_info` object, which is just a wrapper around a list with
//' fields `Q` and `pi_tcag`. The former has the rate matrix, and the latter
//' has the equilibrium nucleotide densities for "T", "C", "A", and "G", respectively.
//' Access the rate matrix for a `sub_model_info` object named `x` via `x$Q` and
//' densities via `x$pi_tcag`.
//'
//' @examples
//' # Same substitution rate for all types:
//' Q_JC69 <- sub_JC69(lambda = 0.1)
//'
//' # Transitions 2x more likely than transversions:
//' Q_K80 <- sub_K80(alpha = 0.2, beta = 0.1)
//'
//' # Same as above, but incorporating equilibrium frequencies
//' sub_HKY85(pi_tcag = c(0.1, 0.2, 0.3, 0.4),
//'           alpha = 0.2, beta = 0.1)
//'
NULL_ENTRY;






//' @describeIn sub_models TN93 model.
//'
//' @param pi_tcag Vector of length 4 indicating the equilibrium distributions of
//'     T, C, A, and G respectively. Values must be >= 0, and
//'     they are forced to sum to 1.
//' @param alpha_1 Substitution rate for T <-> C transition.
//' @param alpha_2 Substitution rate for A <-> G transition.
//' @param beta Substitution rate for transversions.
//'
//' @export
//'
//[[Rcpp::export]]
List sub_TN93(std::vector<double> pi_tcag,
              const double& alpha_1,
              const double& alpha_2,
              const double& beta) {

    // Check that pi_tcag is proper size, no values < 0, at least one value > 0.
    vec_check(pi_tcag, "pi_tcag", true, 4);

    if (alpha_1 < 0) str_stop({"\nFor the TN93 model, `alpha_1` should be >= 0."});
    if (alpha_2 < 0) str_stop({"\nFor the TN93 model, `alpha_2` should be >= 0."});
    if (beta < 0) str_stop({"\nFor the TN93 model, `beta` should be >= 0."});

    // Standardize pi_tcag first:
    double pi_sum = std::accumulate(pi_tcag.begin(), pi_tcag.end(), 0.0);
    for (double& d : pi_tcag) d /= pi_sum;

    arma::mat Q(4,4);
    Q.fill(beta);
    Q.submat(arma::span(0,1), arma::span(0,1)).fill(alpha_1);
    Q.submat(arma::span(2,3), arma::span(2,3)).fill(alpha_2);
    for (uint64 i = 0; i < 4; i++) Q.col(i) *= pi_tcag[i];

    // Reset diagonals to zero
    Q.diag().fill(0.0);

    List out = List::create(_["Q"] = Q, _["pi_tcag"] = pi_tcag);

    out.attr("class") = "sub_model_info";

    return out;
}



//' @describeIn sub_models JC69 model.
//'
//' @param lambda Substitution rate for all possible substitutions.
//'
//' @export
//'
//'
//[[Rcpp::export]]
List sub_JC69(double lambda) {

    if (lambda < 0) str_stop({"\nFor the JC69 model, `lambda` should be >= 0."});

    std::vector<double> pi_tcag(4, 0.25);
    lambda *= 4; // bc it's being multiplied by pi_tcag

    List out = sub_TN93(pi_tcag, lambda, lambda, lambda);

    return out;
}


//' @describeIn sub_models K80 model.
//'
//' @param alpha Substitution rate for transitions.
//' @inheritParams sub_TN93
//'
//' @export
//'
//[[Rcpp::export]]
List sub_K80(double alpha,
             double beta) {

    if (alpha < 0) str_stop({"\nFor the K80 model, `alpha` should be >= 0."});
    if (beta < 0) str_stop({"\nFor the K80 model, `beta` should be >= 0."});

    std::vector<double> pi_tcag(4, 0.25);
    alpha *= 4;  // bc they're being multiplied by pi_tcag
    beta *= 4;  // bc they're being multiplied by pi_tcag

    List out = sub_TN93(pi_tcag, alpha, alpha, beta);

    return out;
}


//' @describeIn sub_models F81 model.
//'
//' @inheritParams sub_TN93
//'
//' @export
//'
//[[Rcpp::export]]
List sub_F81(const std::vector<double>& pi_tcag) {

    List out = sub_TN93(pi_tcag, 1, 1, 1);

    return out;
}


//' @describeIn sub_models HKY85 model.
//'
//'
//' @inheritParams sub_TN93
//' @inheritParams sub_K80
//'
//' @export
//'
//[[Rcpp::export]]
List sub_HKY85(const std::vector<double>& pi_tcag,
               const double& alpha,
               const double& beta) {

    if (alpha < 0) str_stop({"\nFor the HKY85 model, `alpha` should be >= 0."});
    if (beta < 0) str_stop({"\nFor the HKY85 model, `beta` should be >= 0."});

    List out = sub_TN93(pi_tcag, alpha, alpha, beta);

    return out;
}


//' @describeIn sub_models F84 model.
//'
//'
//' @inheritParams sub_TN93
//' @inheritParams sub_K80
//' @param kappa The transition/transversion rate ratio.
//'
//' @export
//'
//[[Rcpp::export]]
List sub_F84(const std::vector<double>& pi_tcag,
             const double& beta,
             const double& kappa) {

    if (beta < 0) str_stop({"\nFor the F84 model, `beta` should be >= 0."});
    if (kappa < 0) str_stop({"\nFor the F84 model, `kappa` should be >= 0."});

    double pi_y = pi_tcag[0] + pi_tcag[1];
    double pi_r = pi_tcag[2] + pi_tcag[3];

    double alpha_1 = (1 + kappa / pi_y) * beta;
    double alpha_2 = (1 + kappa / pi_r) * beta;

    List out = sub_TN93(pi_tcag, alpha_1, alpha_2, beta);

    return out;
}




//' @describeIn sub_models GTR model.
//'
//' @inheritParams sub_TN93
//' @param abcdef A vector of length 6 that contains the off-diagonal elements
//'     for the substitution rate matrix.
//'     See `vignette("sub-models")` for how the values are ordered in the matrix.
//'
//' @export
//'
//[[Rcpp::export]]
List sub_GTR(std::vector<double> pi_tcag,
             const std::vector<double>& abcdef) {

    // Check that pi_tcag is proper size, no values < 0, at least one value > 0.
    vec_check(pi_tcag, "pi_tcag", true, 4);
    // Check that abcdef is proper size, no values < 0.
    vec_check(abcdef, "abcdef", false, 6);

    // Standardize pi_tcag first:
    double pi_sum = std::accumulate(pi_tcag.begin(), pi_tcag.end(), 0.0);
    for (double& d : pi_tcag) d /= pi_sum;

    arma::mat Q(4, 4, arma::fill::zeros);

    // Filling in non-diagonals
    uint64 k = 0;
    for (uint64 i = 0; i < 3; i++) {
        for (uint64 j = i+1; j < 4; j++) {
            Q(i,j) = abcdef[k];
            Q(j,i) = abcdef[k];
            k++;
        }
    }
    for (uint64 i = 0; i < 4; i++) Q.col(i) *= pi_tcag[i];

    List out = List::create(_["Q"] = Q, _["pi_tcag"] = pi_tcag);

    out.attr("class") = "sub_model_info";

    return out;

}




//' @describeIn sub_models UNREST model.
//'
//'
//' @param Q Matrix of substitution rates for "T", "C", "A", and "G", respectively.
//'     Item `Q[i,j]` is the rate of substitution from nucleotide `i` to nucleotide `j`.
//'     Do not include indel rates here!
//'     Values on the diagonal are calculated inside the function so are ignored.
//'
//' @export
//'
//'
//[[Rcpp::export]]
List sub_UNREST(arma::mat Q) {

    /*
     This function also fills in a vector of equilibrium frequencies for each nucleotide.
     This calculation has to be done for this model only because it uses separate
     values for each non-diagonal cell and doesn't use equilibrium frequencies for
     creating the matrix.
     It does this by solving for πQ = 0 by finding the left eigenvector of Q that
     corresponds to the eigenvalue closest to zero.
     */

    if (Q.n_rows != 4 || Q.n_cols != 4) {
        str_stop({"\nIn UNREST model, the matrix `Q` should have 4 rows and 4 columns."});
    }

    /*
     Standardize `Q` matrix
     */
    // Make sure diagonals are set to zero so summing by row works
    Q.diag().fill(0.0);
    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums *= -1;
    Q.diag() = rowsums;

    /*
     Estimate pi_tcag:
     */
    std::vector<double> pi_tcag(4);
    arma::cx_vec eigvals;
    arma::cx_mat eigvecs;

    arma::eig_gen(eigvals, eigvecs, Q.t());

    arma::vec vals = arma::abs(arma::real(eigvals));
    arma::mat vecs = arma::real(eigvecs);

    uint64 i = arma::as_scalar(arma::find(vals == arma::min(vals), 1));

    arma::vec left_vec = vecs.col(i);
    double sumlv = arma::accu(left_vec);

    for (uint64 i = 0; i < 4; i++) pi_tcag[i] = left_vec(i) / sumlv;

    /*
     Assemble final output:
    */
    // Reset diagonal to zero for later steps
    Q.diag().fill(0.0);

    List out = List::create(_["Q"] = Q, _["pi_tcag"] = pi_tcag);

    out.attr("class") = "sub_model_info";

    return out;
}

