
/*
 ****************************
 ****************************

 Methods for the `site_var` R function

 ****************************
 ****************************
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp>   // pcg prng
#include <vector>               // vector class
#include <string>               // string class
#include <random>               // gamma_distribution

#include "jackalope_types.h"       // integer types
#include "pcg.h"                // pcg seeding



using namespace Rcpp;






//' Fill matrix of Gamma-region end points and Gamma values.
//'
//' @param gamma_mat The gamma matrix to fill.
//' @param gammas_x_sizes The value of `sum(gamma[i] * region_size[i])` to fill in.
//'     This value is used to later determine (in fxn `make_gamma_mats`) the
//'     mean gamma value across the whole genome, which is then used to make sure that
//'     the overall mean is 1.
//' @param seq_size_ Length of the focal sequence.
//' @param region_size_ Size of each Gamma region.
//' @param shape The shape parameter for the Gamma distribution from which
//'     Gamma values will be derived.
//' @param invariant Proportion of invariant regions.
//' @param eng A random number generator.
//'
//'
//' @noRd
//'
void fill_gamma_mat_(arma::mat& gamma_mat,
                     double& gammas_x_sizes,
                     const uint64& seq_size_,
                     const uint64& region_size_,
                     const double& shape,
                     const double& invariant,
                     pcg64& eng) {

    // Number of gamma values needed:
    uint64 n_gammas = static_cast<uint64>(std::ceil(
        static_cast<double>(seq_size_) / static_cast<double>(region_size_)));

    // Initialize output matrix
    gamma_mat.set_size(n_gammas, 2);

    // Number of regions that are invariant
    uint64 n_inv = std::round(n_gammas * invariant);
    // Sample indices for which ones are invariant:
    std::deque<uint64> iv_inds;
    if (n_inv > 0) {
        iv_inds = as<std::deque<uint64>>(
            Rcpp::sample(n_gammas, n_inv, false, R_NilValue, false));
        std::sort(iv_inds.begin(), iv_inds.end());
    }


    if (shape > 0) {

        // Initialize Gamma distribution
        std::gamma_distribution<double> distr(shape, shape);

        /*
        Fill matrix.
        Note that I'm setting `start_` to 1 bc the SequenceGammas constructor assumes
        1-based indexing.
        I'm doing it this way to make it more straightforward if someone wants to pass
        their own matrix directly from R (since R obviously uses 1-based).
        */
        double gamma_, size_;
        uint64 end_;
        for (uint64 i = 0, start_ = 1; i < n_gammas; i++, start_ += region_size_) {

            if (!iv_inds.empty() && (i == iv_inds.front())) {
                gamma_ = 0;
                iv_inds.pop_front();
            } else gamma_ = distr(eng);

            end_ = start_ + region_size_ - 1;
            if (i == n_gammas - 1) end_ = seq_size_;

            gamma_mat(i,0) = end_;
            gamma_mat(i,1) = gamma_;

            size_ = end_ - start_ + 1;
            gammas_x_sizes += (size_ * gamma_);

        }

    } else {

        double gamma_, size_;
        uint64 end_;
        for (uint64 i = 0, start_ = 1; i < n_gammas; i++, start_ += region_size_) {

            if (!iv_inds.empty() && (i == iv_inds.front())) {
                gamma_ = 0;
                iv_inds.pop_front();
            } else gamma_ = 1;

            end_ = start_ + region_size_ - 1;
            if (i == n_gammas - 1) end_ = seq_size_;

            gamma_mat(i,0) = end_;
            gamma_mat(i,1) = gamma_;

            size_ = end_ - start_ + 1;
            gammas_x_sizes += (size_ * gamma_);

        }

    }



    return;
}


//' Make matrices of Gamma-region end points and Gamma values for multiple sequences.
//'
//' @param seq_sizes Lengths of the sequences in the genome.
//' @param region_size_ Size of each Gamma region.
//' @param shape The shape parameter for the Gamma distribution from which
//'     Gamma values will be derived.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::field<arma::mat> make_gamma_mats(const std::vector<uint64>& seq_sizes,
                                       const uint64& region_size_,
                                       const double& shape,
                                       const double& invariant) {

    if (region_size_ == 0) {
        stop("\nIn internal function make_gamma_mats, region_size_ == 0");
    }

    pcg64 eng = seeded_pcg();

    uint64 n_seqs = seq_sizes.size();

    /*
     These are used later to make sure the average gamma value across the whole
     genome is 1
     */
    // Total genome size:
    double total_size = static_cast<double>(
        std::accumulate(seq_sizes.begin(), seq_sizes.end(), 0));
    // Sum of each gamma value times the region size it represents:
    double gammas_x_sizes = 0;

    arma::field<arma::mat> gamma_mats(n_seqs);
    for (uint64 i = 0; i < n_seqs; i++) {
        fill_gamma_mat_(gamma_mats(i), gammas_x_sizes, seq_sizes[i],
                        region_size_, shape, invariant, eng);
    }

    /*
     Making sure mean value of gammas across genome equals exactly 1.
     (Only necessary if shape is > 0.)
     */
    if (shape > 0) {
        double mean_gamma = gammas_x_sizes / total_size;
        for (uint64 i = 0; i < n_seqs; i++) gamma_mats(i).col(1) /= mean_gamma;
    }

    return gamma_mats;

}

//' Check input Gamma matrices for proper # columns and end points.
//'
//' @param mats List of matrices to check.
//' @param seq_sizes Vector of sequences sizes for all sequences.
//'
//' @return A length-2 vector of potential error codes and the index (1-based indexing)
//'     to which matrix was a problem.
//'
//' @noRd
//'
//[[Rcpp::export]]
void check_gamma_mats(const std::vector<arma::mat>& mats,
                      const std::vector<uint64>& seq_sizes) {

    std::string err_msg = "\nIf providing custom matrices for the ";
    err_msg += "`mats` argument to the `site_var` function, ";
    err_msg += "all matrices ";

    bool error = false;

    for (uint64 i = 0; i < mats.size(); i++) {

        const arma::mat& gamma_mat(mats[i]);

        // There needs to be a column of end points and a column of gamma values
        if (gamma_mat.n_cols != 2) {
            err_msg += "need to have 2 columns, one for end positions, one for gamma ";
            err_msg += "distances.";
            error = true;
            break;
        }
        // Since the input matrix should be 1-based indexing, make sure there are no 0s:
        if (arma::any(gamma_mat.col(0) <= 0)) {
            err_msg += "should only have values > 0 in the first column, ";
            err_msg += "which is where the end points should be.";
            error = true;
            break;
        }
        // Sampling weights < 0 makes no sense
        if (arma::any(gamma_mat.col(1) < 0)) {
            err_msg += "should only have values >= 0 in the second column, ";
            err_msg += "which is where the mutation-rate weights should be.";
            error = true;
            break;
        }
        // Make sure there are no repeat ends points
        arma::vec unq_ends = arma::unique(gamma_mat.col(0));
        if (unq_ends.n_elem != gamma_mat.n_rows) {
            err_msg += "should contain no duplicate end points (in the first column).";
            error = true;
            break;
        }
        // Non-integer numbers in the end points doesn't make sense:
        arma::vec trunc_ends = arma::trunc(gamma_mat.col(0));
        if (arma::any(gamma_mat.col(0) != trunc_ends)) {
            err_msg += "should contain only whole numbers as end points ";
            err_msg += "(i.e., in the first column).";
            error = true;
            break;
        }
        // The last end point should be the end of the sequence:
        uint64 last_end = static_cast<uint64>(gamma_mat.col(0).max());
        if (last_end != seq_sizes[i]) {
            err_msg += "need to have a maximum end point (in the first column) ";
            err_msg += "equal to the size of the associated sequence.";
            error = true;
            break;
        }
    }

    if (error) throw(Rcpp::exception(err_msg.c_str(), false));

    return;
}
