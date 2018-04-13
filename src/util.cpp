
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <random> // C++11 sampling distributions
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include <RcppArmadilloExtensions/sample.h>

#include "gemino_types.h"
#include "sequence_classes.h" // new classes
#include "alias.h" // alias sampling
#include "util.h"

using namespace Rcpp;


/*
 ========================================================================================
 ========================================================================================

 Standard sampling

 ========================================================================================
 ========================================================================================
 */




// Returns an approximately normal distribution, except that it's discrete and cannot be
// be less than `min_len`. Used in `rando_seqs` and `create_genome` below.

arma::ivec non_zero_norm(uint N, double mean_len, double sd_len, double min_len) {

    arma::vec vec_a = arma::randn(N);
    vec_a = (vec_a * sd_len) + mean_len;
    vec_a = arma::round(vec_a);
    vec_a.elem(arma::find(vec_a < min_len)).fill(min_len);

    arma::ivec vec_b = arma::conv_to<arma::ivec>::from(vec_a);

    return vec_b;
}







/*
 ========================================================================================
 ========================================================================================

 Random sequences

 ========================================================================================
 ========================================================================================
 */

// For one char vector, no weighting

std::string cpp_rando_seq(uint len) {
    std::string bases = "ACGT";
    std::string out_str(len, 'x');
    for (uint i = 0; i < len; i++) {
        uint ind = static_cast<uint>(R::unif_rand() * 4);
        out_str[i] = bases[ind];
    }
    return out_str;
}


//' Create random sequences.
//'
//' Function to create random sequences into character vector format or for a new
//' reference genome object.
//'
//' Note that this function will not return sequences smaller than 10bp.
//'
//' @param N Number of sequences.
//' @param mean_len Mean length of each sequence.
//' @param mean_len Standard deviation of lengths.
//' @param min_len Minimum length of any sequence. Defaults to \code{1}.
//'
//' @return Character vector of sequence strings.
//'
//' @export
//'
//' @examples
//' randos <- rando_seqs(1000)
//'
//[[Rcpp::export]]
std::vector<std::string> rando_seqs(int N,
                                    double mean_len = 100,
                                    double sd_len = 0,
                                    double min_len = 1) {

    arma::ivec n_vec = non_zero_norm(N, mean_len, sd_len, min_len);

    std::string bases = "ACGT";

    uint base_index;
    std::vector<std::string> out_vec(N, "");

    for (int i = 0; i < N; i++) {
        int& n_i = n_vec[i];
        for (int j = 0; j < n_i; j++) {
            base_index = static_cast<uint>(R::unif_rand() * 4);
            out_vec[i] += bases[base_index];
        }
        Rcpp::checkUserInterrupt();
    }

    return out_vec;

}



//' Create RefGenome pointer based on nucleotide equilibrium frequencies.
//'
//' Function to create random sequences into character vector format for or for a new
//' reference genome object.
//'
//' Note that this function will never return empty sequences.
//'
//' @param N Number of sequences.
//' @param mean_ Mean for the gamma distribution for sequence sizes.
//' @param sd_ Standard deviation for the gamma distribution for sequence sizes.
//' @param Vector of nucleotide equilibrium frequencies for "A", "C", "G", and "T",
//'     respectively. Defaults to \code{rep(0.25, 4)}.
//'
//' @return Character vector of sequence strings.
//'
//' @export
//'
//' @examples
//'
//' genome <- create_genome(10, 100e6, 10e6, pis = c(0.1, 0.2, 0.3, 0.4))
//'
//[[Rcpp::export]]
SEXP create_genome(const uint& N,
                   const double& mean_,
                   const double& sd_,
                   NumericVector equil_freqs = NumericVector(0),
                   uint n_cores = 1) {

    if (sd_ <= 0) stop("sd_ must be > 0, otherwise the function `create_genome` hangs.");

    if (equil_freqs.size() == 0) equil_freqs = NumericVector(4, 0.25);

    // Converting to STL format
    std::vector<double> pis = as<std::vector<double>>(equil_freqs);

    // Generate seeds for random number generators (1 per core)
    const std::vector<uint> seeds = as<std::vector<uint>>(
        Rcpp::runif(n_cores, 0, static_cast<double>(sitmo::prng_engine::max())));

    // Alias-sampling object
    const alias_FL fl(pis);

    // Creating output object
    XPtr<RefGenome> ref_xptr(new RefGenome(N), true);
    // Reference to be used for manipulating
    RefGenome& ref(*ref_xptr);

    // parameters for creating the gamma distribution
    const double gamma_shape = (mean_ * mean_) / (sd_ * sd_);
    const double gamma_scale = (sd_ * sd_) / mean_;


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
        std::gamma_distribution<double> distr(gamma_shape, gamma_scale);

        // Parallelize the Loop
        #ifdef _OPENMP
        #pragma omp for schedule(static)
        #endif
        for (uint i = 0; i < N; i++) {
            RefSequence& seq(ref.sequences[i]);
            seq.name = "seq" + std::to_string(i);
            // Get length of output sequence:
            uint len = static_cast<uint>(distr(engine));
            if (len < 1) len = 1;
            // Sample sequence:
            seq.nucleos = std::string(len, 'x');
            alias_sample_str(seq.nucleos, fl, engine);
        }

    #ifdef _OPENMP
    }
    #endif

    for (uint i = 0; i < N; i++) {
        ref.total_size += ref[i].size();
    }

    return ref_xptr;
}





/*
 ========================================================================================
 ========================================================================================

 Sequence info

 ========================================================================================
 ========================================================================================
 */


//' GC proportion of a single string.
//'
//'
//' @param sequence String for a single sequence.
//'
//' @return Proportion of sequence that's a `'G'` or `'C'`.
//'
//' @noRd
//'
double gc_prop(const std::string& sequence) {
    int total_seq = sequence.size();
    double total_gc = 0;
    for (int i = 0; i < total_seq; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
// ... overloaded for portion of a string
double gc_prop(const std::string& sequence, const uint& start, const uint& stop) {
    int total_seq = stop - start + 1;
    double total_gc = 0;
    for (int i = start; i <= stop; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            total_gc += 1;
        }
    }
    double gc_prop = total_gc / total_seq;
    return gc_prop;
}
