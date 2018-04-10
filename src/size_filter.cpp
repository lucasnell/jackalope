#include <RcppArmadillo.h>
#include <vector>  // Vector class
#include <cmath>  // std::exp, std::log
#include <random>  // uniform_real_distribution
#include <sitmo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
#include <progress_bar.hpp>


using namespace Rcpp;

#include "gemino_types.h"  // integer types, SequenceSet, VariantSet



// Convert digestion cut site locations into fragment sizes
std::vector<uint> digests_to_sizes(const std::vector<uint>& one_digests, const uint& total_size) {

    uint n = one_digests.size();

    std::vector<uint> out(n+1);
    out[0] = one_digests[0];

    for (uint i = 1; i < n; i++) {
        out[i] = one_digests[i] - one_digests[(i-1)];
    }

    out[n] = total_size - one_digests[(n-1)];

    return out;
}



// Density for lognormal distribution
std::vector<double> cpp_dlnorm(const std::vector<uint>& frag_sizes, double mu, double sigma) {
    std::vector<double> out_vec(frag_sizes.size());
    double x, out, e_exp;
    for (uint i = 0; i < frag_sizes.size(); i++) {
        x = static_cast<double>(frag_sizes[i]);
        out = 1 / (x * sigma * std::sqrt(2 * M_PI));
        e_exp = std::log(x) - mu;
        e_exp *= e_exp;
        e_exp /= (2 * std::exp(2 * std::log(sigma)));
        e_exp *= -1;
        out *= std::exp(e_exp);
        out_vec[i] = out;
    }
    return out_vec;
}
// Overloading for a single fragment size
double cpp_dlnorm(const uint& frag_size, double mu, double sigma) {
    double x, out, e_exp;
    x = static_cast<double>(frag_size);
    out = 1 / (x * sigma * std::sqrt(2 * M_PI));
    e_exp = std::log(x) - mu;
    e_exp *= e_exp;
    e_exp /= (2 * std::exp(2 * std::log(sigma)));
    e_exp *= -1;
    out *= std::exp(e_exp);
    return out;
}



// Gets frequencies of fragment sizes, and combines all sizes > 1,000 bp into one
// category
// For a bin with index i, it represents fragment sizes in the range...
// [(i*increment)+1, ((i*increment))+increment]
// Except for the last item, which is all fragments > 1,000 bp

// Note: Fragment sizes == 0 will be skipped entirely

std::vector<double> freq_vector(const std::vector<uint>& frag_sizes, uint increment = 1) {

    if ((1000 % increment) != 0) stop("Increment should be a factor of 1,000.");
    uint n = 1000 / increment;
    n++;  // Adding extra item for fragments > 1,000 bp

    std::vector<double> out(n);
    uint ind;
    for (uint i = 0; i < frag_sizes.size(); i++) {
        if (frag_sizes[i] == 0) continue;
        ind = (frag_sizes[i] - 1) / increment;
        if (ind > n-1) ind = n-1;
        out[ind] += 1;
    }
    return out;
}




// Do all the filtering for one scaffold's fragment sizes

std::vector<uint> size_filter_scaff(const std::vector<uint>& frag_sizes,
                               sitmo::prng_engine& eng,
                               const uint& increment) {


    std::vector<uint> output(0);
    uint n_frags = frag_sizes.size();
    if (n_frags == 0) return output;

    // Table of frequencies by fragment size
    std::vector<double> f_tab = freq_vector(frag_sizes, increment);
    // This is just casting the n_frags into a double
    double N = n_frags;

    double mean_rho = std::exp(0.005101 - 1.004 * std::log(increment));

    // pr_seq is proportional to P(Sequenced)
    // pr_avail is proportional to P(Available to be sequenced)
    // prob is P(Sequenced | Available)
    // U holds value for a uniform random number in range [0,1)
    double pr_seq, pr_avail, prob, U;
    // For drawing U values
    std::uniform_real_distribution<double> distr(0, 1);
    // Index for the focal fragment size
    uint frag_ind;

    for (uint i = 0; i < n_frags; i++) {
        pr_seq = cpp_dlnorm(frag_sizes[i], 4.724, 0.6950);
        frag_ind = (frag_sizes[i] - 1) / increment;
        if (frag_ind > (f_tab.size() - 1)) frag_ind = f_tab.size() - 1;
        pr_avail = f_tab[frag_ind] / N;
        if (pr_avail == 0) Rcout << "pr_aval == 0" << std::endl;
        prob = pr_seq / pr_avail;
        prob *= (0.3119 / mean_rho);
        U = distr(eng);
        if (U < prob) output.push_back(i); // Keeping in C++ indices
    }

    return output;
}





// ====================================================================================
// ====================================================================================
// ====================================================================================

//      Size filter reference

// ====================================================================================
// ====================================================================================
// ====================================================================================



// Weighted sampling of a vector of fragment sizes for one reference's scaffold.
//
// The weighting is based on the combined probability densities from the predicted and
// sequenced fragments in maize from Elshire et al. (2011;
// doi:10.1371/journal.pone.0019379).
//
// See size_filter.Rmd vignette for more info
//
// @param frag_sizes A vector of fragment sizes.
// @param eng A sitmo pseudo-random number generator object.
// @param increment The width of fragment sizes to group together when calculating
//      the availability probability.
//
// @return A vector of indices for the fragment sizes that were selected.
//
std::vector<uint> size_filter_ref_scaff(const uint& scaff_index,
                                   const std::vector<std::vector<uint>>& all_digests,
                                   const XPtr<SequenceSet>& reference,
                                   sitmo::prng_engine& eng,
                                   uint increment) {


    std::vector<uint> output(0);
    const std::vector<uint>& one_digests = all_digests[scaff_index];
    if (one_digests.size() == 0) return output;
    const uint& scaff_size = reference->seq_sizes[scaff_index];

    std::vector<uint> frag_sizes = digests_to_sizes(one_digests, scaff_size);

    // Adjusting merged scaffold's first fragment length
    if (reference->merged && scaff_index > 0) {
        if (all_digests[(scaff_index-1)].size() > 0) {
            frag_sizes[0] += (reference->seq_sizes[(scaff_index-1)] -
                all_digests[(scaff_index-1)].back());
        } else {
            uint ind = scaff_index;
            while (ind > 0) {
                if (all_digests[(ind-1)].size() == 0) break;
                frag_sizes[0] += reference->seq_sizes[(ind-1)];
                ind--;
            }
            if (ind > 0) {
                frag_sizes[0] += (reference->seq_sizes[(ind-1)] -
                    all_digests[(ind-1)].back());
            }
        }
    }

    // uint n_frags = frag_sizes.size();

    output = size_filter_scaff(frag_sizes, eng, increment);

    return output;
}



//' Weighted sampling of a digested \code{dna_set}'s fragment sizes.
//'
//' See \code{size_filter.Rmd} vignette for more information on the sampling scheme and
//' the justification for it.
//'
//' Each index present in the output from this function indicates to keep the
//' fragment \emph{preceding} the cut site at that index.
//' For example, if you have the following conditions:
//' \enumerate{
//'   \item The second scaffold for your focal genome consists of the following
//'     sequence: \code{"AACCGGTT"}.
//'   \item The digestion cut sites for this scaffold are at positions \code{0},
//'     \code{3}, and \code{5}.
//'   \item The output from this function for this scaffold is \code{1}.
//' }
//'
//' ... then you'd keep the fragment from index \code{0} to index \code{2}: \code{"AAC"}.
//'
//' \emph{Note:} All indices are in C++ 0-based format.
//'
//'
//' @param all_digests A list of vectors of digestion cut sites, one vector per scaffold.
//' @param reference A \code{dna_set} object's \code{sequence_set} field, which
//'     represents genome information.
//' @param n_cores Number of cores to use for processing. This argument is ignored if
//'     OpenMP is not enabled. Defaults to 1.
//' @param increment The width of fragment sizes to group together when calculating
//'      the availability probability. Defaults to 1.
//'
//' @return A list of vectors of indices for the fragment sizes that were selected.
//'
//[[Rcpp::export]]
std::vector<std::vector<uint>> filter_reference_frags(const std::vector<std::vector<uint>>& all_digests,
                                            const XPtr<SequenceSet>& reference,
                                            const uint& n_cores = 1,
                                            const uint& increment = 1) {

    const std::vector<uint> seeds = as<std::vector<uint>>(runif(n_cores, 0, 4294967296));

    const uint n_scaffs = reference->sequences.size();

    std::vector<std::vector<uint>> out(n_scaffs, std::vector<uint>(0));

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seed = seeds[omp_get_thread_num()];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint i = 0; i < n_scaffs; i++) {
        out[i] = size_filter_ref_scaff(i, all_digests, reference, engine, increment);
    }
    #ifdef _OPENMP
    }
    #endif

    return out;
}





// ====================================================================================
// ====================================================================================
// ====================================================================================

//      Size filter variants

// ====================================================================================
// ====================================================================================
// ====================================================================================







std::vector<uint> size_filter_var_scaff(
        const uint& scaff_index,
        const uint& var_index,
        const std::vector< std::vector< std::vector< uint > > >& all_digests,
        const XPtr<SequenceSet>& reference,
        const XPtr<VariantSet>& variant_set,
        sitmo::prng_engine& eng,
        uint increment) {

    // This variant's digests
    const std::vector<std::vector<uint>>& var_digests = all_digests[var_index];

    std::vector<uint> output(0);
    const std::vector<uint>& one_digests = var_digests[scaff_index];
    if (one_digests.size() == 0) return output;
    const OneVariant& var_inf = variant_set->variant_info[var_index];
    const uint& scaff_size = var_inf.scaffold_lengths[scaff_index];

    std::vector<uint> frag_sizes = digests_to_sizes(one_digests, scaff_size);

    // Adjusting merged scaffold's first fragment length
    if (reference->merged && scaff_index > 0) {
        // var_inf.scaffold_lengths[scaff_index];
        if (var_digests[(scaff_index-1)].size() > 0) {
            frag_sizes[0] += (var_inf.scaffold_lengths[(scaff_index-1)] -
                var_digests[(scaff_index-1)].back());
        } else {
            uint ind = scaff_index;
            while (ind > 0) {
                if (var_digests[(ind-1)].size() == 0) break;
                frag_sizes[0] += var_inf.scaffold_lengths[(ind-1)];
                ind--;
            }
            if (ind > 0) {
                frag_sizes[0] += (var_inf.scaffold_lengths[(ind-1)] -
                    var_digests[(ind-1)].back());
            }
        }
    }

    // uint n_frags = frag_sizes.size();

    output = size_filter_scaff(frag_sizes, eng, increment);

    return output;
}




//' Weighted sampling of a digested \code{variants}'s fragment sizes.
//'
//' See \code{size_filter.Rmd} vignette for more information on the sampling scheme and
//' the justification for it.
//'
//' Each index present in the output from this function indicates to keep the
//' fragment \emph{preceding} the cut site at that index.
//' For example, if you have the following conditions:
//' \enumerate{
//'   \item The second scaffold for your focal genome consists of the following
//'     sequence: \code{"AACCGGTT"}.
//'   \item The digestion cut sites for this scaffold are at positions \code{0},
//'     \code{3}, and \code{5}.
//'   \item The output from this function for this scaffold is \code{1}.
//' }
//'
//' ... then you'd keep the fragment from index \code{0} to index \code{2}: \code{"AAC"}.
//'
//' \emph{Note:} All indices are in C++ 0-based format.
//'
//'
//' @param all_digests A list of lists of vectors of digestion cut sites, one
//'     vector for each unique scaffold--variant combination.
//' @param reference A \code{variants} object's \code{reference} field, which
//'     represents genome information.
//' @param variant_set A \code{variants} object's \code{variant_set} field, which
//'     represents information on variants' deviations from the reference genome.
//' @param n_cores Number of cores to use for processing. This argument is ignored if
//'     OpenMP is not enabled. Defaults to 1.
//' @param increment The width of fragment sizes to group together when calculating
//'      the availability probability. Defaults to 1.
//'
//' @return A list of lists of vectors of indices for the fragment sizes that
//'     were selected.
//'
//[[Rcpp::export]]
std::vector< std::vector< std::vector< uint > > > filter_variants_frags(
        const std::vector< std::vector< std::vector< uint > > >& all_digests,
        const XPtr<SequenceSet>& reference,
        const XPtr<VariantSet>& variant_set,
        const uint n_cores = 1,
        const uint increment = 1) {

    const std::vector<uint> seeds = as<std::vector<uint>>(runif(n_cores, 0, 4294967296));

    const uint n_scaffs = reference->sequences.size();
    const uint n_vars = variant_set->variant_info.size();

    std::vector<std::vector<std::vector<uint>>> out(
            n_vars, std::vector<std::vector<uint>>(
                    n_scaffs, std::vector<uint>(0)));


    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    uint active_seed;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    active_seed = seeds[omp_get_thread_num()];
    #else
    active_seed = seeds[0];
    #endif

    sitmo::prng_engine engine(active_seed);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint j = 0; j < n_vars; j++) {
        for (uint i = 0; i < n_scaffs; i++) {
            out[j][i] = size_filter_var_scaff(i, j, all_digests, reference,
                                              variant_set, engine, increment);
        }
    }
    #ifdef _OPENMP
    }
    #endif

    return out;
}



