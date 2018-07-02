
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
#include "sequence_classes.h"  // Var* and Ref* classes
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

    XPtr<std::vector<std::vector<PhyloTree>>> trees(
            new std::vector<std::vector<PhyloTree>>(n_seqs)
    );

    for (uint32 i = 0; i < n_seqs; i++) {
        const List& seq_phylo_info(genome_phylo_info[i]);
        std::vector<PhyloTree>& trees_i((*trees)[i]);
        uint32 n_trees = seq_phylo_info.size();
        trees_i = std::vector<PhyloTree>(n_trees);
        for (uint32 j = 0; j < n_trees; j++) {
            const List& phylo_info(seq_phylo_info[j]);
            std::vector<double> branch_lens = as<std::vector<double>>(
                phylo_info["branch_lens"]);
            arma::Mat<uint32> edges = as<arma::Mat<uint32>>(phylo_info["edges"]);
            std::vector<std::string> labels = as<std::vector<std::string>>(
                phylo_info["labels"]);
            uint32 start = as<uint32>(phylo_info["start"]);
            sint64 end = as<sint64>(phylo_info["end"]);
            trees_i[j] = PhyloTree(branch_lens, edges, labels, start, end);
        }
    }

    return trees;
}








