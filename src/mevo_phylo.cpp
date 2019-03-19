
/*
 ********************************************************

 Methods for molecular evolution using phylogenies

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <progress.hpp>  // for the progress bar


#include "jackal_types.h"  // integer types
#include "seq_classes_var.h"  // Var* classes
#include "mevo.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "mevo_phylo.h"


using namespace Rcpp;






/*
`T` should be `MutationSampler` or `ChunkMutationSampler`.
*/
template <typename T>
XPtr<VarSet> PhyloInfo<T>::evolve_seqs(
        SEXP& ref_genome_ptr,
        SEXP& sampler_base_ptr,
        const std::vector<arma::mat>& gamma_mats,
        uint32 n_cores,
        const bool& show_progress) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    XPtr<T> sampler_base(sampler_base_ptr);

    // Extract tip labels from the first tree:
    const PhyloTree& first_tree(phylo_one_seqs[0].trees[0]);
    std::vector<std::string> var_names = first_tree.tip_labels;

    XPtr<VarSet> var_set(new VarSet(*ref_genome, var_names), true);

#ifndef _OPENMP
    n_cores = 1;
#endif

    uint32 n_seqs = var_set->reference->size();
    uint64 total_seq = ref_genome->total_size;

    Progress prog_bar(total_seq, show_progress);
    std::vector<int> status_codes(n_cores, 0);

    if (n_seqs != gamma_mats.size()) {
        std::string err_msg = "\ngamma_mats must be of same length as # sequences in ";
        err_msg += "reference";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }
    if (n_seqs != phylo_one_seqs.size()) {
        std::string err_msg = "\n# tips in phylo. info must be of same length as ";
        err_msg += "# sequences in reference genome";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }

    for (uint32 i = 0; i < n_seqs; i++) {
        if (gamma_mats[i](gamma_mats[i].n_rows-1,0) != (*var_set)[0][i].size()) {
            std::string err_msg = "\nGamma matrices must have max values equal to ";
            err_msg += "the respective sequence's length.\n";
            err_msg += "This error occurred on Gamma matrix number ";
            err_msg += std::to_string(i+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }
    }


    // Generate seeds for random number generators (1 RNG per core)
    const std::vector<std::vector<uint64>> seeds = mc_seeds(n_cores);

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_cores) if (n_cores > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
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

        if (status_code != 0) continue;

        PhyloOneSeq<T>& seq_phylo(phylo_one_seqs[i]);

        const arma::mat& gamma_mat(gamma_mats[i]);

        // Set values for variant info and sampler:
        seq_phylo.set_samp_var_info(*var_set, *sampler_base, i, gamma_mat);

        // Evolve the sequence using the seq_phylo object:
        status_code = seq_phylo.evolve(eng, prog_bar);

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
        }
    }

    return var_set;

}







//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees(const List& genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<PhyloInfo<MutationSampler>> all_seqs_xptr(
            new PhyloInfo<MutationSampler>(genome_phylo_info)
    );

    return all_seqs_xptr;
}


//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' Same as above, but chunked.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees_chunk(const List& genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<PhyloInfo<ChunkMutationSampler>> all_seqs_xptr(
            new PhyloInfo<ChunkMutationSampler>(genome_phylo_info)
    );

    return all_seqs_xptr;
}





//' Evolve all sequences in a reference genome.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP evolve_seqs(
        SEXP& ref_genome_ptr,
        SEXP& sampler_base_ptr,
        SEXP& phylo_info_ptr,
        const std::vector<arma::mat>& gamma_mats,
        const uint32& n_cores,
        const bool& show_progress) {

    XPtr<PhyloInfo<MutationSampler>> phylo_info(phylo_info_ptr);

    XPtr<VarSet> var_set = phylo_info->evolve_seqs(
        ref_genome_ptr, sampler_base_ptr,
        gamma_mats, n_cores, show_progress);

    return var_set;
}

//' Same as above, but using chunks.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP evolve_seqs_chunk(
        SEXP& ref_genome_ptr,
        SEXP& sampler_base_ptr,
        SEXP& phylo_info_ptr,
        const std::vector<arma::mat>& gamma_mats,
        const uint32& n_cores,
        const bool& show_progress) {

    XPtr<PhyloInfo<ChunkMutationSampler>> phylo_info(phylo_info_ptr);

    XPtr<VarSet> var_set = phylo_info->evolve_seqs(
        ref_genome_ptr, sampler_base_ptr,
        gamma_mats, n_cores, show_progress);

    return var_set;
}






