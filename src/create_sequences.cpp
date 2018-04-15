
/*
 ********************************************************

 Functions to create new sequences.

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <random> // C++11 sampling distributions
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "gemino_types.h"  // integer types
#include "sequence_classes.h" // RefGenome, RefSequence classes
#include "alias.h" // alias sampling

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
OuterClass create_sequences_(const uint& n_seqs,
                             const double& len_mean,
                             const double& len_sd,
                             NumericVector equil_freqs,
                             const uint& n_cores) {

    if (equil_freqs.size() == 0) equil_freqs = NumericVector(4, 0.25);

    // Converting to STL format
    std::vector<double> pis = as<std::vector<double>>(equil_freqs);

    // Generate seeds for random number generators (1 per core)
    const std::vector<uint> seeds = as<std::vector<uint>>(
        Rcpp::runif(n_cores, 0, static_cast<double>(sitmo::prng_engine::max())));

    // Alias-sampling object
    const alias_FL fl(pis);

    // Creating output object
    OuterClass seqs_out(n_seqs);

    // parameters for creating the gamma distribution
    const double gamma_shape = (len_mean * len_mean) / (len_sd * len_sd);
    const double gamma_scale = (len_sd * len_sd) / len_mean;


    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint active_thread = omp_get_thread_num();
    active_seed = seeds[active_thread];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);
    // Gamma distribution to be used for size selection (doi: 10.1093/molbev/msr011):
    std::gamma_distribution<double> distr;
    if (len_sd > 0) {
        distr = std::gamma_distribution<double>(gamma_shape, gamma_scale);
    }

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint i = 0; i < n_seqs; i++) {
        InnerClass& seq(seqs_out[i]);

        // Get length of output sequence:
        uint len;
        if (len_sd > 0) {
            len = static_cast<uint>(distr(engine));
            if (len < 1) len = 1;
        } else len = len_mean;
        // Sample sequence:
        seq.resize(len, 'x');
        alias_sample_str<InnerClass>(seq, fl, engine);
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
//'     If set to `<= 0`, all sequences will be the same length. Defaults to `0`.
//' @param equil_freqs Vector of nucleotide equilibrium frequencies for
//'     "A", "C", "G", and "T", respectively. Defaults to `rep(0.25, 4)`.
//' @param n_cores Number of cores to use via OpenMP.
//'
//'
//' @return External pointer to a `RefGenome` C++ object.
//'
//' @export
//'
//' @examples
//'
//' genome <- create_genome(10, 100e6, 10e6, equil_freqs = c(0.1, 0.2, 0.3, 0.4))
//'
//[[Rcpp::export]]
SEXP create_genome(const uint& n_seqs,
                   const double& len_mean,
                   const double& len_sd = 0,
                   NumericVector equil_freqs = NumericVector(0),
                   const uint& n_cores = 1) {

    RefGenome ref;

    ref = create_sequences_<RefGenome, RefSequence>(
        n_seqs, len_mean, len_sd, equil_freqs, n_cores);

    for (uint i = 0; i < n_seqs; i++) {
        ref.total_size += ref[i].size();
        ref[i].name = "seq" + std::to_string(i);
    }

    XPtr<RefGenome> ref_xptr(&ref);

    return ref_xptr;
}




//' `rando_seqs` creates random sequences as a character vector.
//'
//'
//' @inheritParams create_genome
//'
//' @return Character vector of sequence strings.
//'
//' @describeIn create_genome
//'
//' @export
//'
//' @examples
//' randos <- rando_seqs(10, 1000, 10)
//'
//[[Rcpp::export]]
std::vector<std::string> rando_seqs(const uint& n_seqs,
                                    const double& len_mean,
                                    const double& len_sd = 0,
                                    NumericVector equil_freqs = NumericVector(0),
                                    const uint& n_cores = 1) {

    std::vector<std::string> ref = create_sequences_<std::vector<std::string>,
                std::string>(n_seqs, len_mean, len_sd, equil_freqs, n_cores);

    return ref;
}

