//
// Alter reference genome sequences
//



#include <RcppArmadillo.h>

#include <algorithm>  // random_shuffle
#include <deque>
#include <string>
#include <vector>
#include <progress.hpp>  // for the progress bar

#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "seq_classes_ref.h"  // Ref* classes
#include "jackalope_types.h"  // integer types
#include "alias_sampler.h"  // alias string sampler
#include "util.h"  // clear_memory, thread_check


using namespace Rcpp;


namespace alter_scaffs {
    // wrapper around R's RNG such that we get a uniform distribution over
    // [0,n) as required by the STL algorithm
    // (see http://gallery.rcpp.org/articles/stl-random-shuffle/)
    inline int rand_wrapper(const int n) { return std::floor(unif_rand()*n); }
}



// ======================================================================================
// ======================================================================================

//  Merge sequences

// ======================================================================================
// ======================================================================================



//' Merge a reference genome into a single sequence.
//'
//'
//' @param ref_genome_ptr An external pointer (R class \code{externalptr}) to a
//'     \code{RefGenome} class in C++ (the full class in C++ is
//'     \code{Rcpp::XPtr<RefGenome>}).
//'
//' @return Nothing. Changes are made in place.
//'
//' @name merge_sequences
//'
//' @noRd
//'
//[[Rcpp::export]]
void merge_sequences_cpp(SEXP ref_genome_ptr) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefSequence>& seqs(ref_genome->sequences);

    // Shuffling ref_genome info.
    std::random_shuffle(seqs.begin(), seqs.end(), alter_scaffs::rand_wrapper);

    // Merging the back sequences to the first one:
    std::string& nts(seqs.front().nucleos);
    ref_genome->old_names.push_back(seqs.front().name);
    seqs.front().name = "MERGE";
    uint64 i = seqs.size() - 1;
    while (seqs.size() > 1) {
        nts += seqs[i].nucleos;
        ref_genome->old_names.push_back(seqs[i].name);
        --i;
        seqs.pop_back();
    }
    // clear memory in string
    clear_memory<std::string>(nts);
    // clear memory in deque
    clear_memory<std::deque<RefSequence>>(seqs);

    ref_genome->merged = true;

    return;
}









// ======================================================================================
// ======================================================================================

//  Filter sequences

// ======================================================================================
// ======================================================================================

//' Filter reference genome sequences by size or for a proportion of total nucleotides.
//'
//'
//' @inheritParams ref_genome_ptr merge_sequences
//' @param min_seq_size Integer minimum sequence size to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//' @param out_seq_prop Numeric proportion of total sequence to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//'
//' @return Nothing. Changes are made in place.
//'
//' @name filter_sequences
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void filter_sequences_cpp(SEXP ref_genome_ptr,
                          const uint64& min_seq_size = 0,
                          const double& out_seq_prop = 0) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefSequence>& seqs(ref_genome->sequences);

    // Checking for sensible inputs
    if (out_seq_prop <= 0 && min_seq_size == 0) {
        stop("Specify > 0 for min_seq_size or out_seq_prop");
    }
    if (out_seq_prop > 0 && min_seq_size > 0) {
        stop("Specify > 0 for min_seq_size OR out_seq_prop");
    }
    if (out_seq_prop > 1) stop("out_seq_prop must be between 0 and 1");

    // Sorting sequence set by size (largest first)
    std::sort(seqs.begin(), seqs.end(), std::greater<RefSequence>());

    // Index that will point to the first sequence to be deleted
    uint64 i = 0;
    // Keeping track of total genome size after filtering
    double out_seq = 0;

    if (min_seq_size > 0) {
        if (seqs.back().size() >= min_seq_size) return;
        if (seqs[i].size() < min_seq_size) {
            str_stop({"Desired minimum scaffold size is too large. None found. ",
                     "The minimum size is ", std::to_string(seqs[i].size())});
        }
        // after below, `iter` points to the first sequence smaller than the minimum
        while (seqs[i].size() >= min_seq_size) {
            out_seq += static_cast<double>(seqs[i].size());
            ++i;
        }
    } else {
        // Changing total_size to double so I don't have to worry about integer division
        // being a problem
        double total_seq = static_cast<double>(ref_genome->total_size);
        out_seq = static_cast<double>(seqs[i].size());
        while (out_seq / total_seq < out_seq_prop) {
            ++i;
            out_seq += static_cast<double>(seqs[i].size());
        }
        // Getting `i` to point to the first item to be deleted:
        ++i;
    }

    // Erasing using `iter`
    if (i < seqs.size()) {
        seqs.erase(seqs.begin() + i, seqs.end());
        // clear memory:
        clear_memory<std::deque<RefSequence>>(seqs);
    }

    ref_genome->total_size = static_cast<uint64>(out_seq);

    return;
}





// ======================================================================================
// ======================================================================================

//  Replace Ns with random sequences

// ======================================================================================
// ======================================================================================

//' Replace Ns with randome nucleotides.
//'
//'
//' @return Nothing. Changes are made in place.
//'
//' @name replace_Ns_cpp
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void replace_Ns_cpp(SEXP ref_genome_ptr,
                    const std::vector<double>& pi_tcag,
                    uint64 n_threads,
                    const bool& show_progress) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    // Check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

    const uint64 n_seqs = ref_genome->size();

    // Progress bar
    Progress prog_bar(ref_genome->total_size, show_progress);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads) if (n_threads > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint64 active_thread = omp_get_thread_num();
#else
    uint64 active_thread = 0;
#endif
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    // Samples for nucleotides:
    AliasStringSampler<std::string> sampler("TCAG", pi_tcag);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint64 i = 0; i < n_seqs; i++) {
        if (prog_bar.is_aborted() || prog_bar.check_abort()) continue;
        RefSequence& seq(ref_genome->sequences[i]);
        for (char& c : seq.nucleos) {
            if (c == 'N') c = sampler.sample(eng);
        }
        prog_bar.increment(seq.size());
    }

#ifdef _OPENMP
}
#endif

    return;
}
