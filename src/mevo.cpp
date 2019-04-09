
#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <progress.hpp>  // for the progress bar
#ifdef _OPENMP
#include <omp.h>  // omp
#endif


#include "jackalope_types.h"
#include "mevo.h"
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class

using namespace Rcpp;





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





//' Used below to directly make a MutationTypeSampler
//'
//' @noRd
//'
MutationTypeSampler make_type_sampler(const arma::mat& Q,
                                      const std::vector<double>& pi_tcag,
                                      const std::vector<double>& insertion_rates,
                                      const std::vector<double>& deletion_rates) {

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

    // Type and insertion samplers:
    MutationTypeSampler type(probs, mut_lengths);

    return type;
}




//' Add mutations at segregating sites for one sequence from coalescent simulation output.
//'
//' @noRd
//'
void add_one_seq_sites(VarSet& var_set,
                       const RefGenome& ref_genome,
                       const uint32& seq_i,
                       const arma::mat& ss_i,
                       MutationTypeSampler& type_sampler,
                       AliasStringSampler<std::string>& insert_sampler,
                       pcg64& eng) {

    uint32 pos;
    std::string nts; // <-- for insertions
    /*
     Going from back so we don't have to update positions constantly, and so we can
     use the char from the reference genome directly.
     */
    for (uint32 k = 0; k < ss_i.n_rows; k++) {
        uint32 i = ss_i.n_rows - 1 - k;
        pos = ss_i(i, 0);
        MutationInfo mut = type_sampler.sample(ref_genome[seq_i][pos], eng);
        if (mut.length == 0) {
            for (uint32 j = 1; j < ss_i.n_cols; j++) {
                if (ss_i(i,j) == 1) {
                    var_set[j-1][seq_i].add_substitution(mut.nucleo, pos);
                }
            }
        } else if (mut.length > 0) {
            nts.resize(mut.length);  // resize nts on insertion len
            insert_sampler.sample(nts, eng);  // fill w/ random nucleotides
            for (uint32 j = 1; j < ss_i.n_cols; j++) {
                if (ss_i(i,j) == 1) {
                    var_set[j-1][seq_i].add_insertion(nts, pos);
                }
            }
        } else {
            sint64 pos_ = static_cast<sint64>(pos);
            sint64 size_ = static_cast<sint64>(var_set.min_size(seq_i));
            if (pos_ - mut.length > size_) {
                mut.length = static_cast<sint32>(pos_-size_);
            }
            uint32 del_size = std::abs(mut.length);
            for (uint32 j = 1; j < ss_i.n_cols; j++) {
                if (ss_i(i,j) == 1) {
                    var_set[j-1][seq_i].add_deletion(del_size, pos);
                }
            }
        }
    }

    return;
}





//' Add mutations at segregating sites from coalescent simulation output.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP add_coal_sites_cpp(SEXP& ref_genome_ptr,
                        const std::vector<std::string>& var_names,
                        const std::vector<arma::mat>& seg_sites,
                        const arma::mat& Q,
                        const std::vector<double>& pi_tcag,
                        const std::vector<double>& insertion_rates,
                        const std::vector<double>& deletion_rates,
                        uint32 n_threads,
                        const bool& show_progress) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    // Initialize new VarSet object
    XPtr<VarSet> var_set(new VarSet(*ref_genome, var_names), true);

#ifndef _OPENMP
    n_threads = 1;
#endif

    const uint32 n_seqs = ref_genome->size();
    const uint64 total_seq = ref_genome->total_size;

    Progress prog_bar(total_seq, show_progress);
    std::vector<int> status_codes(n_threads, 0);

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
{
#endif

    // Type and insertion samplers:
    MutationTypeSampler type = make_type_sampler(Q, pi_tcag, insertion_rates,
                                                 deletion_rates);
    AliasStringSampler<std::string> insert(mevo::bases, pi_tcag);

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
#else
    uint32 active_thread = 0;
#endif
    int& status_code(status_codes[active_thread]);
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint32 i = 0; i < n_seqs; i++) {

        if (prog_bar.is_aborted() || prog_bar.check_abort()) status_code = -1;
        if (status_code != 0) continue;

        add_one_seq_sites(*var_set, *ref_genome, i, seg_sites[i], type, insert, eng);

        prog_bar.increment((*ref_genome)[i].size());

    }

#ifdef _OPENMP
}
#endif

    for (const int& status_code : status_codes) {
        if (status_code == -1) {
            std::string warn_msg = "\nThe user interrupted phylogenetic evolution. ";
            warn_msg += "Note that changes occur in place, so your variants have ";
            warn_msg += "already been partially added.";
            Rcpp::warning(warn_msg.c_str());
            break;
        }
    }

    return var_set;

}
