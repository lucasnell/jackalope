
/*
 ********************************************************

 Functions to create new sequences.

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <random> // C++11 sampling distributions
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "alias_sampler.h" // alias sampling
#include "pcg.h" // pcg::max, mt_seeds, seeded_pcg
#include "util.h" // thread_check

using namespace Rcpp;









/*
 ========================================================================================
 ========================================================================================

 Random sequences

 ========================================================================================
 ========================================================================================
 */



/*
 Template that does most of the work for creating new sequences for the following two
 functions when `len_sd > 0`.
 Classes `OuterClass` and `InnerClass` can be `std::vector<std::string>` and
 `std::string` or `RefGenome` and `RefSequence`.
 No other combinations are guaranteed to work.
 */

template <typename OuterClass, typename InnerClass>
OuterClass create_sequences_(const uint64& n_seqs,
                             const double& len_mean,
                             const double& len_sd,
                             const std::vector<double>& pi_tcag,
                             uint64 n_threads) {

    // Check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

    Progress prog_bar(n_seqs, false); // just use as way to check for abort

    // Alias-sampling object
    const AliasSampler sampler(pi_tcag);

    // Creating output object
    OuterClass seqs_out(n_seqs);

    // parameters for creating the gamma distribution
    const double gamma_shape = (len_mean * len_mean) / (len_sd * len_sd);
    const double gamma_scale = (len_sd * len_sd) / len_mean;


    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
    #ifdef _OPENMP
    uint64 active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
    #else
    active_seeds = seeds[0];
    #endif

    pcg64 engine = seeded_pcg(active_seeds);
    // Gamma distribution to be used for size selection (doi: 10.1093/molbev/msr011):
    std::gamma_distribution<double> distr;
    if (len_sd > 0) {
        distr = std::gamma_distribution<double>(gamma_shape, gamma_scale);
    }

    std::string bases_ = alias_sampler::bases;

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint64 i = 0; i < n_seqs; i++) {

        if (prog_bar.is_aborted() || prog_bar.check_abort()) continue;

        InnerClass& seq(seqs_out[i]);

        // Get length of output sequence:
        uint64 len;
        if (len_sd > 0) {
            len = static_cast<uint64>(distr(engine));
            if (len < 1) len = 1;
        } else len = len_mean;
        // Sample sequence:
        seq.reserve(len);
        for (uint64 j = 0; j < len; j++) {
            uint64 k = sampler.sample(engine);
            seq.push_back(bases_[k]);
        }
    }

    #ifdef _OPENMP
    }
    #endif

    return seqs_out;
}





//' Create `RefGenome` pointer based on nucleotide equilibrium frequencies.
//'
//' Function to create random sequences for a new reference genome object.
//'
//' Note that this function will never return empty sequences.
//'
//' @param n_seqs Number of sequences.
//' @param len_mean Mean for the gamma distribution for sequence sizes.
//' @param len_sd Standard deviation for the gamma distribution for sequence sizes.
//'     If set to `<= 0`, all sequences will be the same length.
//' @param pi_tcag Vector of nucleotide equilibrium frequencies for
//'     "T", "C", "A", and "G", respectively.
//' @param n_threads Number of threads to use via OpenMP.
//'
//'
//' @return External pointer to a `RefGenome` C++ object.
//'
//' @noRd
//'
//' @examples
//'
//'
//[[Rcpp::export]]
SEXP create_genome_cpp(const uint64& n_seqs,
                       const double& len_mean,
                       const double& len_sd,
                       std::vector<double> pi_tcag,
                       const uint64& n_threads) {

    XPtr<RefGenome> ref_xptr(new RefGenome(), true);
    RefGenome& ref(*ref_xptr);

    ref = create_sequences_<RefGenome, RefSequence>(
        n_seqs, len_mean, len_sd, pi_tcag, n_threads);

    for (uint64 i = 0; i < n_seqs; i++) {
        ref.total_size += ref[i].size();
        ref[i].name = "seq" + std::to_string(i);
    }

    return ref_xptr;
}




//' Create random sequences as a character vector.
//'
//' This function is used internally for testing.
//'
//'
//' @inheritParams create_genome
//'
//' @return Character vector of sequence strings.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> rando_seqs(const uint64& n_seqs,
                                    const double& len_mean,
                                    const double& len_sd = 0,
                                    NumericVector pi_tcag = NumericVector(0),
                                    const uint64& n_threads = 1) {

    std::vector<double> pi_tcag_ = as<std::vector<double>>(pi_tcag);
    if (pi_tcag_.size() == 0) pi_tcag_ = std::vector<double>(4, 0.25);

    std::vector<std::string> ref = create_sequences_<std::vector<std::string>,
                std::string>(n_seqs, len_mean, len_sd, pi_tcag_, n_threads);

    return ref;
}

