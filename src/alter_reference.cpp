//
// Alter reference genome chromosomes
//


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>

#include <algorithm>  // random_shuffle, sort
#include <deque>
#include <string>
#include <vector>
#include <progress.hpp>  // for the progress bar

#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "ref_classes.h"  // Ref* classes
#include "jackalope_types.h"  // integer types
#include "alias_sampler.h"  // alias string sampler
#include "util.h"  // clear_memory, thread_check, jlp_shuffle


using namespace Rcpp;




// ======================================================================================
// ======================================================================================

//  Merge chromosomes

// ======================================================================================
// ======================================================================================


/*
 Merge all chromosomes from a reference genome into a single chromosome.
*/
//[[Rcpp::export]]
void merge_all_chromosomes_cpp(SEXP ref_genome_ptr) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefChrom>& chroms(ref_genome->chromosomes);

    pcg64 eng = seeded_pcg();

    // Shuffling ref_genome info.
    jlp_shuffle<std::deque<RefChrom>>(chroms, eng);

    // Merging the back chromosomes to the first one:
    std::string& nts(chroms.front().nucleos);
    ref_genome->old_names.push_back(chroms.front().name);
    chroms.front().name = "MERGE";
    uint64 i = chroms.size() - 1;
    while (chroms.size() > 1) {
        nts += chroms[i].nucleos;
        ref_genome->old_names.push_back(chroms[i].name);
        --i;
        chroms.pop_back();
    }
    // clear memory in string
    clear_memory<std::string>(nts);
    // clear memory in deque
    clear_memory<std::deque<RefChrom>>(chroms);

    ref_genome->merged = true;

    return;
}

/*
 Merge one or more reference chromosomes into a single chromosome.
*/
//[[Rcpp::export]]
void merge_chromosomes_cpp(SEXP ref_genome_ptr,
                           std::deque<uint64> chrom_inds) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefChrom>& chroms(ref_genome->chromosomes);

    // Merging the back chromosomes to the first one:
    RefChrom& chrom(chroms[chrom_inds.front()]);
    std::string& nts(chrom.nucleos);

    for (uint64 i = 1; i < chrom_inds.size(); i++) {
        RefChrom& chrom_i(chroms[chrom_inds[i]]);
        chrom.name += "__";
        chrom.name += chrom_i.name;
        std::string& nts_i(chrom_i.nucleos);
        nts += nts_i;
        // clear memory in string
        nts_i.clear();
        clear_memory<std::string>(nts_i);
    }
    // Go back and remove RefChrom objects:
    chrom_inds.pop_front(); // don't want to remove first one w all the sequence!
    std::sort(chrom_inds.begin(), chrom_inds.end());
    for (auto iter = chrom_inds.rbegin(); iter != chrom_inds.rend(); ++iter) {
        chroms.erase(chroms.begin() + *iter);
    }
    // clear memory in deque
    clear_memory<std::deque<RefChrom>>(chroms);

    return;
}









// ======================================================================================
// ======================================================================================

//  Filter chromosomes

// ======================================================================================
// ======================================================================================

//' Filter reference genome chromosomes by size or for a proportion of total nucleotides.
//'
//'
//' @inheritParams ref_genome_ptr merge_chromosomes
//' @param min_chrom_size Integer minimum chromosome size to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//' @param out_chrom_prop Numeric proportion of total chromosome to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//'
//' @return Nothing. Changes are made in place.
//'
//' @name filter_chromosomes
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void filter_chromosomes_cpp(SEXP ref_genome_ptr,
                          const uint64& min_chrom_size = 0,
                          const double& out_chrom_prop = 0) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefChrom>& chroms(ref_genome->chromosomes);

    // Checking for sensible inputs
    if (out_chrom_prop <= 0 && min_chrom_size == 0) {
        stop("Specify > 0 for min_chrom_size or out_chrom_prop");
    }
    if (out_chrom_prop > 0 && min_chrom_size > 0) {
        stop("Specify > 0 for min_chrom_size OR out_chrom_prop");
    }
    if (out_chrom_prop > 1) stop("out_chrom_prop must be between 0 and 1");

    // Sorting chromosome set by size (largest first)
    std::sort(chroms.begin(), chroms.end(), std::greater<RefChrom>());

    // Index that will point to the first chromosome to be deleted
    uint64 i = 0;
    // Keeping track of total genome size after filtering
    double out_chrom = 0;

    if (min_chrom_size > 0) {
        if (chroms.back().size() >= min_chrom_size) return;
        if (chroms[i].size() < min_chrom_size) {
            str_stop({"Desired minimum chromosome size is too large. None found. ",
                     "The largest chromosome is ", std::to_string(chroms[i].size())});
        }
        // after below, `iter` points to the first chromosome smaller than the minimum
        while (chroms[i].size() >= min_chrom_size) {
            out_chrom += static_cast<double>(chroms[i].size());
            ++i;
        }
    } else {
        // Changing total_size to double so I don't have to worry about integer division
        // being a problem
        double total_chrom = static_cast<double>(ref_genome->total_size);
        out_chrom = static_cast<double>(chroms[i].size());
        while (out_chrom / total_chrom < out_chrom_prop) {
            ++i;
            out_chrom += static_cast<double>(chroms[i].size());
        }
        // Getting `i` to point to the first item to be deleted:
        ++i;
    }

    // Erasing using `iter`
    if (i < chroms.size()) {
        chroms.erase(chroms.begin() + i, chroms.end());
        // clear memory:
        clear_memory<std::deque<RefChrom>>(chroms);
    }

    ref_genome->total_size = static_cast<uint64>(out_chrom);

    return;
}





// ======================================================================================
// ======================================================================================

//  Replace Ns with random chromosomes

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

    const uint64 n_chroms = ref_genome->size();

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
    for (uint64 i = 0; i < n_chroms; i++) {
        if (prog_bar.is_aborted() || prog_bar.check_abort()) continue;
        RefChrom& chrom(ref_genome->chromosomes[i]);
        for (char& c : chrom.nucleos) {
            if (c == 'N') c = sampler.sample(eng);
        }
        prog_bar.increment(chrom.size());
    }

#ifdef _OPENMP
}
#endif

    return;
}
