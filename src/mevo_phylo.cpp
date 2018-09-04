
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


#include "gemino_types.h"  // integer types
#include "seq_classes_var.h"  // Var* classes
#include "mevo.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "mevo_rate_matrices.h"
#include "mevo_phylo.h"


using namespace Rcpp;






//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees(List genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<std::vector<PhyloOneSeq<MutationSampler>>> all_seqs_xptr(
            new std::vector<PhyloOneSeq<MutationSampler>>(n_seqs)
    );
    std::vector<PhyloOneSeq<MutationSampler>>& all_seqs(*all_seqs_xptr);

    for (uint32 i = 0; i < n_seqs; i++) {
        fill_one_seq_<MutationSampler>(genome_phylo_info, i, all_seqs);
    }

    return all_seqs_xptr;
}


//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' Same as above, but chunked.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees_chunk(List genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<std::vector<PhyloOneSeq<ChunkMutationSampler>>> all_seqs_xptr(
            new std::vector<PhyloOneSeq<ChunkMutationSampler>>(n_seqs)
    );
    std::vector<PhyloOneSeq<ChunkMutationSampler>>& all_seqs(*all_seqs_xptr);

    for (uint32 i = 0; i < n_seqs; i++) {
        fill_one_seq_<ChunkMutationSampler>(genome_phylo_info, i, all_seqs);
    }

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
        const std::vector<uint32>& seq_inds,
        const std::vector<arma::mat>& gamma_mats,
        const bool& show_progress) {

    XPtr<VarSet> var_set =
        evolve_seqs_<MutationSampler>(ref_genome_ptr,
                                      sampler_base_ptr,
                                      phylo_info_ptr,
                                      seq_inds, gamma_mats, show_progress);

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
        const std::vector<uint32>& seq_inds,
        const std::vector<arma::mat>& gamma_mats,
        const bool& show_progress) {

    XPtr<VarSet> var_set =
        evolve_seqs_<ChunkMutationSampler>(ref_genome_ptr,
                                           sampler_base_ptr,
                                           phylo_info_ptr,
                                           seq_inds, gamma_mats, show_progress);

    return var_set;
}






