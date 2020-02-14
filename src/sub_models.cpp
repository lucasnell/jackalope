
/*
 ******************************************************************************
 ******************************************************************************

 Create rate matrices for multiple molecular-evolution models:
 TN93 and its special cases (JC69, K80, F81, HKY85, and F84), plus GTR and UNREST

 ******************************************************************************
 ******************************************************************************
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <cmath> // exp, pow, remainder

#include "jackalope_types.h" // integer types
#include "util.h" // str_stop



using namespace Rcpp;







/*
 ======================================================================================
 ======================================================================================

 +Gamma model functions

 ======================================================================================
 ======================================================================================
 */


/*
 Incomplete Gamma function
 */
inline double incG(const double& a, const double& z) {
    return R::pgamma(z, a, 1.0, 0, 0) * R::gammafn(a);
}




/*
 Mean of truncated Gamma distribution
 From http://dx.doi.org/10.12988/astp.2013.310125.
 As in that paper, b > 0 is the scale and c > 0 is the shape.
*/
double trunc_Gamma_mean(const double& b, const double& c,
                        const double& xl, const double& xu) {

    // Ran the following in Mathematica to find out that this part goes to
    // zero as xu goes to infinity:
    // > Limit[b^(-c + 1) Exp[-x/b] x^c, x -> Infinity, Assumptions -> c in real && b > 0]
    // So if xu is Inf, then we set this to zero:
    double k_;
    if (xu == arma::datum::inf) {
        k_ = 0;
    } else {
        k_ = std::exp(-1.0 * xu / b) * std::pow(b, 1-c) * std::pow(xu, c);
    }
    double k = c / (
        b * incG(1+c, xl/b) - b * incG(1+c, xu/b) +
            k_ -
            std::exp(-1.0 * xl / b) * std::pow(b, 1.0 - c) * std::pow(xl, c)
    );
    double z = -(b * b) * k * (- incG(1+c, xl / b) + incG(1+c, xu / b));
    return z;
}

// Create a vector of Gamma values for a discrete Gamma distribution.
void discrete_gamma(std::vector<double>& gammas,
                    const uint32& k,
                    const double& shape) {

    if (shape <= 0 || k <= 1) {
        gammas.push_back(1.0);
        return;
    }

    gammas.reserve(k);

    double scale = 1 / shape;
    double d_k = 1.0 / static_cast<double>(k);

    double p_cutoff = d_k;
    double xl = 0, xu = 0;

    for (uint32 i = 0; i < k; i++) {
        xl = xu;
        if (p_cutoff < 1) {
            xu = R::qgamma(p_cutoff, shape, scale, 1, 0);
        } else xu = arma::datum::inf;
        gammas.push_back(trunc_Gamma_mean(scale, shape, xl, xu));
        p_cutoff += d_k;
    }

    return;

}


// Info to calculate P(t) for TN93 model and its special cases
void Pt_info(std::vector<double> pi_tcag,
             const double& alpha_1,
             const double& alpha_2,
             const double& beta,
             arma::mat& U,
             arma::mat& Ui,
             arma::vec& L) {


    const double& pi_t(pi_tcag[0]);
    const double& pi_c(pi_tcag[1]);
    const double& pi_a(pi_tcag[2]);
    const double& pi_g(pi_tcag[3]);

    double pi_y = pi_t + pi_c;
    double pi_r = pi_a + pi_g;

    U = {
        {1,     1 / pi_y,       0,              pi_c / pi_y},
        {1,     1 / pi_y,       0,              -pi_t / pi_y},
        {1,     -1 / pi_r,      pi_g / pi_r,    0},
        {1,     -1 / pi_r,      -pi_a / pi_r,   0}};

    Ui = {
        {pi_t,          pi_c,           pi_a,           pi_g},
        {pi_t * pi_r,   pi_c * pi_r,    -pi_a * pi_y,   -pi_g * pi_y},
        {0,             0,              1,              -1},
        {1,             -1,             0,              0}};

    L = {0, -beta, -(pi_r * alpha_2 + pi_y * beta), -(pi_y * alpha_1 + pi_r * beta)};

    return;
}

// Info to calculate P(t) for GTR model
void Pt_info(const arma::mat& Q,
             arma::mat& U,
             arma::mat& Ui,
             arma::vec& L) {

    // Special case when all Q are zeros:
    if (arma::all(arma::vectorise(Q) == 0)) {

        U = {{1,    2,      0.0,    0.5},
            {1,     2,      0.0,    -0.5},
            {1,     -2,     0.5,    0.0},
            {1,     -2,     -0.5,   0.0}};
        Ui = {{0.250,     0.250,      0.250,      0.250},
            {0.125,     0.125,      -0.125,     -0.125},
            {0.000,     0.000,      1.000,      -1.000},
            {1.000,     -1.000,     0.000,       0.000}};

        L = arma::vec(4, arma::fill::zeros);

        return;
    }


    arma::cx_vec L_;
    arma::cx_mat U_;

    arma::eig_gen(L_, U_, Q);

    arma::uvec I = arma::sort_index(arma::abs(arma::real(L_)), "descend");
    L_ = L_(I);
    U_ = U_.cols(I);

    L = arma::real(L_);
    U = arma::real(U_);

    Ui = U.i();

    return;

}














/*
 ======================================================================================
 ======================================================================================

 Substitution model functions

 ======================================================================================
 ======================================================================================
 */



// Scale rate matrix so that total rate is `mu`
void scale_Q(arma::mat& Q,
             const std::vector<double>& pi_tcag,
             const double& mu) {
    if (mu > 0) {
        arma::vec q(pi_tcag);
        q %= Q.diag();
        double mu_now = std::abs(arma::accu(q));
        if (mu_now != 0) Q *= (mu / mu_now);
    }
    return;
}






//[[Rcpp::export]]
List sub_TN93_cpp(const double& mu,
                  std::vector<double> pi_tcag,
                  const double& alpha_1,
                  const double& alpha_2,
                  const double& beta,
                  const double& gamma_shape,
                  const uint32& gamma_k,
                  const double& invariant) {

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
    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums *= -1;
    Q.diag() = rowsums;

    // Scaling if desired:
    scale_Q(Q, pi_tcag, mu);

    // Now getting vector of Gammas (which is the vector { 1 } if gamma_shape <= 0)
    std::vector<double> gammas;
    discrete_gamma(gammas, gamma_k, gamma_shape);

    // Extract info for P(t)
    arma::field<arma::mat> Q_vec(gammas.size());
    arma::field<arma::mat> U(gammas.size());
    arma::field<arma::mat> Ui(gammas.size());
    arma::field<arma::vec> L(gammas.size());
    for (uint32 i = 0; i < gammas.size(); i++) {
        Q_vec(i) = Q * gammas[i];
        Pt_info(pi_tcag, alpha_1 * gammas[i], alpha_2 * gammas[i], beta * gammas[i],
                U(i), Ui(i), L(i));
    }

    List out = List::create(_["Q"] = Q_vec,
                            _["pi_tcag"] = pi_tcag,
                            _["U"] = U,
                            _["Ui"] = Ui,
                            _["L"] = L,
                            _["gammas"] = gammas,
                            _["invariant"] = invariant,
                            _["model"] = "TN93");


    return out;
}





//[[Rcpp::export]]
List sub_GTR_cpp(const double& mu,
                 std::vector<double> pi_tcag,
                 const std::vector<double>& abcdef,
                 const double& gamma_shape,
                 const uint32& gamma_k,
                 const double& invariant) {

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

    // Filling in diagonals
    arma::vec rowsums = arma::sum(Q, 1);
    rowsums *= -1;
    Q.diag() = rowsums;

    // Scaling if desired:
    scale_Q(Q, pi_tcag, mu);

    // Now getting vector of Gammas (which is the vector { 1 } if gamma_shape <= 0)
    std::vector<double> gammas;
    discrete_gamma(gammas, gamma_k, gamma_shape);

    // Extract info for P(t)
    arma::field<arma::mat> Q_vec(gammas.size());
    arma::field<arma::mat> U(gammas.size());
    arma::field<arma::mat> Ui(gammas.size());
    arma::field<arma::vec> L(gammas.size());
    for (uint32 i = 0; i < gammas.size(); i++) {
        Q_vec(i) = Q * gammas[i];
        Pt_info(Q_vec(i), U(i), Ui(i), L(i));
    }

    List out = List::create(_["Q"] = Q_vec,
                            _["pi_tcag"] = pi_tcag,
                            _["U"] = U,
                            _["Ui"] = Ui,
                            _["L"] = L,
                            _["gammas"] = gammas,
                            _["invariant"] = invariant,
                            _["model"] = "GTR");

    return out;

}





//[[Rcpp::export]]
List sub_UNREST_cpp(const double& mu,
                    arma::mat Q,
                    const double& gamma_shape,
                    const uint32& gamma_k,
                    const double& invariant) {

    /*
     This function also fills in a vector of equilibrium frequencies for each nucleotide.
     This calculation has to be done for this model only because it uses separate
     values for each non-diagonal cell and doesn't use equilibrium frequencies for
     creating the matrix.
     It does this by solving for πQ = 0 by finding the left eigenvector of Q that
     corresponds to the eigenvalue closest to zero.
     */

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

    // Scaling if desired:
    scale_Q(Q, pi_tcag, mu);

    // Now getting vector of Gammas (which is the vector { 1 } if gamma_shape <= 0)
    std::vector<double> gammas;
    discrete_gamma(gammas, gamma_k, gamma_shape);

    // Info for P(t) calculation
    arma::field<arma::mat> Q_vec(gammas.size());
    arma::field<arma::mat> U;
    arma::field<arma::mat> Ui;
    arma::field<arma::vec> L;
    // Check to see if the eigenvalues are real:
    arma::cx_vec L_;
    arma::cx_mat U_;
    arma::eig_gen(L_, U_, Q);
    bool all_real = arma::all(arma::imag(L_) == 0) &&
        arma::all(arma::vectorise(arma::imag(U_)) == 0);
    // If they're real, fill U, Ui, and L using GTR method.
    // Otherwise they'll be empty, so we'll use repeated matrix squaring.
    if (all_real) {
        U.set_size(gammas.size());
        Ui.set_size(gammas.size());
        L.set_size(gammas.size());
    }
    for (uint32 i = 0; i < gammas.size(); i++) {
        Q_vec(i) = Q * gammas[i];
        if (all_real) {
            Pt_info(Q_vec(i), U(i), Ui(i), L(i));
        }
    }


    List out = List::create(_["Q"] = Q_vec,
                            _["pi_tcag"] = pi_tcag,
                            _["U"] = U,
                            _["Ui"] = Ui,
                            _["L"] = L,
                            _["gammas"] = gammas,
                            _["invariant"] = invariant,
                            _["model"] = "UNREST");

    return out;
}

