
#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class


#include "jackalope_types.h"
#include "mutator_location.h"
#include "mutator_type.h"
#include "mutator.h"

using namespace Rcpp;

// Add mutation and return the change in the sequence rate that results
double MutationSampler::mutate(pcg64& eng) {
    uint32 pos = sample_location(eng, 0, 0, false);
    // uint32 pos = runif_01(eng) * var_seq->size();
    char c = var_seq->get_nt(pos);
    MutationInfo m = sample_type(c, eng);
    double rate_change;
    if (m.length == 0) {
        rate_change = location.substitution_rate_change(m.nucleo, pos);
        var_seq->add_substitution(m.nucleo, pos);
    } else {
        if (m.length > 0) {
            std::string nts = new_nucleos(m.length, eng);
            rate_change = location.insertion_rate_change(nts, pos);
            var_seq->add_insertion(nts, pos);
        } else {
            sint64 pos_ = static_cast<sint64>(pos);
            sint64 size_ = static_cast<sint64>(var_seq->size());
            if (pos_ - m.length > size_) m.length = static_cast<sint32>(pos_-size_);
            uint32 del_size = std::abs(m.length);
            rate_change = location.deletion_rate_change(m.length, pos);
            var_seq->add_deletion(del_size, pos);
        }
        // Update Gamma region bounds:
        location.update_gamma_regions(m.length, pos);

    }
    return rate_change;
}

/*
 Overloaded for only mutating within a range.
 It also updates `end` if an indel occurs in the range.
 Make sure to keep checking for situation where `end < start` (i.e., sequence section
 is empty).
 `// ***` mark difference between this and previous `mutate` versions
 */
double MutationSampler::mutate(pcg64& eng, const uint32& start, sint64& end) {
    if (end < 0) stop("end is negative in MutationSampler.mutate");
    uint32 pos = sample_location(eng, start, static_cast<uint32>(end), true);  // ***
    // uint32 pos = (runif_01(eng) * (static_cast<uint32>(end) - start + 1)) + start;
    char c = var_seq->get_nt(pos);
    MutationInfo m = sample_type(c, eng);
    double rate_change;
    if (m.length == 0) {
        rate_change = location.substitution_rate_change(m.nucleo, pos);
        var_seq->add_substitution(m.nucleo, pos);
    } else {
        if (m.length > 0) {
            std::string nts = new_nucleos(m.length, eng);
            rate_change = location.insertion_rate_change(nts, pos);
            var_seq->add_insertion(nts, pos);
        } else {
            sint64 pos_ = static_cast<sint64>(pos);
            sint64 size_ = end + 1;  // ***
            if (pos_ - m.length > size_) m.length = static_cast<sint32>(pos_-size_);
            uint32 del_size = std::abs(m.length);
            rate_change = location.deletion_rate_change(m.length, pos);
            var_seq->add_deletion(del_size, pos);
        }
        // Update Gamma region bounds:
        location.update_gamma_regions(m.length, pos);
        // Update end point:
        end += static_cast<sint64>(m.length);  // ***
    }
    return rate_change;
}



// Wrapper to make mutation sampler available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const arma::mat& Q,
                                const std::vector<double>& pi_tcag,
                                const std::vector<double>& insertion_rates,
                                const std::vector<double>& deletion_rates,
                                const uint32& chunk_size) {

    std::vector<std::vector<double>> probs;
    std::vector<sint32> mut_lengths;
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

    MutationRates mr(q_tcag);
    out->location = LocationSampler(mr);

    out->location.change_chunk(chunk_size);

    return out;
}


