
#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class


#include "jackalope_types.h"
#include "mutator_location.h"
#include "mutator_type.h"
#include "mutator.h"

using namespace Rcpp;

// Add mutation and return the change in the sequence rate that results
template <class C>
double OneSeqMutationSampler<C>::mutate(pcg64& eng) {
    uint32 pos = sample_location(eng);
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
template <class C>
double OneSeqMutationSampler<C>::mutate(pcg64& eng, const uint32& start, sint64& end) {
    if (end < 0) stop("end is negative in [Chunk]MutationSampler.mutate");
    uint32 pos = sample_location(eng, start, static_cast<uint32>(end), true);  // ***
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


// Explicit template instantiation
template class OneSeqMutationSampler<LocationSampler>;
template class OneSeqMutationSampler<ChunkLocationSampler>;



// Wrapper to make non-chunked version available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const arma::mat& Q,
                                const std::vector<double>& pi_tcag,
                                const std::vector<double>& insertion_rates,
                                const std::vector<double>& deletion_rates) {

    XPtr<MutationSampler> out =
        make_mutation_sampler_base_<MutationSampler,LocationSampler>(
                Q, pi_tcag, insertion_rates, deletion_rates);

    return out;
}

// Same thing, but with chunks

//[[Rcpp::export]]
SEXP make_mutation_sampler_chunk_base(const arma::mat& Q,
                                      const std::vector<double>& pi_tcag,
                                      const std::vector<double>& insertion_rates,
                                      const std::vector<double>& deletion_rates,
                                      const uint32& chunk_size) {

    XPtr<ChunkMutationSampler> out =
        make_mutation_sampler_base_<ChunkMutationSampler,ChunkLocationSampler>(
                Q, pi_tcag, insertion_rates, deletion_rates);

    out->location.change_chunk(chunk_size);

    return out;
}


