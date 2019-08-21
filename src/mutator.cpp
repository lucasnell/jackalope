

#include "jackalope_types.h"

#include <RcppArmadillo.h>
#include <cmath>  // exp
#include <vector>  // vector class
#include <string>  // string class


#include "mutator_location.h"
#include "mutator_type.h"
#include "mutator.h"

using namespace Rcpp;






//' Changing P(t) matrix with new branch lengths or times.
//'
//'
//' Equivalent to `U %*% diag(exp(L * t)) %*% Ui`
//'
//' @noRd
//'
void Pt_calc(const arma::mat& U,
             const arma::mat& Ui,
             const arma::vec& L,
             const double& t,
             arma::mat& Pt) {

    arma::mat diag_L = arma::diagmat(arma::exp(L * t));

    Pt = U * diag_L * Ui;

    return;
}

//' Calculating P(t) using repeated matrix squaring, for UNREST model only.
//'
//' @noRd
//'
void Pt_calc(const arma::mat& Q,
             const uint32& k,
             const double& t,
             arma::mat& Pt) {

    double m = static_cast<double>(1U<<k);

    Pt = arma::eye<arma::mat>(4, 4) + Q * t / m + 0.5 * (Q * t / m) * (Q * t / m);

    for (uint32 i = 0; i < k; i++) Pt = Pt * Pt;

    return;

}






// Does most of the work of mutating for the below methods (all but location sampling)
inline double MutationSampler::mutate__(pcg64& eng, const uint64& pos, sint64& end) {

    char c = var_chrom->get_nt(pos);
    MutationInfo m = type.sample(c, eng);
    double rate_change;
    if (m.length == 0) {
        rate_change = location.substitution_rate_change(m.nucleo, pos);
        var_chrom->add_substitution(m.nucleo, pos);
    } else {
        if (m.length > 0) {
            std::string nts = new_nucleos(m.length, eng);
            rate_change = location.insertion_rate_change(nts, pos);
            var_chrom->add_insertion(nts, pos);
        } else {
            sint64 pos_ = static_cast<sint64>(pos);
            sint64 size_ = end + 1;
            if (pos_ - m.length > size_) m.length = static_cast<sint64>(pos_-size_);
            uint64 del_size = std::abs(m.length);
            rate_change = location.deletion_rate_change(del_size, pos);
            var_chrom->add_deletion(del_size, pos);
        }
        // Update end point:
        end += static_cast<sint64>(m.length);
    }

    // Update regions, rates, and bounds:
    location.update(rate_change, m.length, pos);

    return rate_change;
}


// Add mutation and return the change in the chromosome rate that results
double MutationSampler::mutate(pcg64& eng) {

    uint64 pos = location.sample(eng);

    if (pos >= var_chrom->size()) {
        Rcout << pos << ' ' << var_chrom->size() << std::endl;
        stop("pos returning too large a pos");
    }

    // Dummy end point for use in mutate__
    sint64 end = var_chrom->size() - 1;

    double rate_change = mutate__(eng, pos, end);

    return rate_change;
}

/*
 Overloaded for only mutating within a range.
 It also updates `end` if an indel occurs in the range.
 Make sure to keep checking for situation where `end < start` (i.e., chromosome section
 is empty).
 */
double MutationSampler::mutate(pcg64& eng, const uint64& start, sint64& end) {

    if (end < 0) stop("end is negative in MutationSampler.mutate");
    uint64 pos = location.sample(eng, start, static_cast<uint64>(end));

    double rate_change = mutate__(eng, pos, end);

    return rate_change;
}



// Wrapper to make mutation sampler available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const arma::mat& Q,
                                const std::vector<double>& pi_tcag,
                                const std::vector<double>& insertion_rates,
                                const std::vector<double>& deletion_rates,
                                const uint64& region_size) {

    std::vector<std::vector<double>> probs;
    std::vector<sint64> mut_lengths;
    std::vector<double> q_tcag;
    /*
     (1) Combine substitution, insertion, and deletion rates into a single vector
     (2) Fill the `q_tcag` vector with mutation rates for each nucleotide
     */
    fill_probs_q_tcag(probs, q_tcag, Q, pi_tcag, insertion_rates, deletion_rates);

    // Now filling in mut_lengths vector
    fill_mut_lengths(mut_lengths, insertion_rates, deletion_rates);

    /*
     Now create and fill output pointer to base sampler:
     */
    XPtr<MutationSampler> out(new MutationSampler());

    out->type = MutationTypeSampler(probs, mut_lengths);
    out->insert = AliasStringSampler<std::string>("TCAG", pi_tcag);

    out->location = LocationSampler(q_tcag, region_size);

    return out;
}


