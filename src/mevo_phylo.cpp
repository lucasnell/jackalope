
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
#include "mevo_phylo.h"


using namespace Rcpp;






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






