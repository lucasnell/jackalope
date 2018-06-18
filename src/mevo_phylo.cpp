
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
#include "molecular_evolution.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "mevo_rate_matrices.h"
#include "mevo_phylo.h"


using namespace Rcpp;



/*
 This function produces `spp_order` in `test_phylo`.
 This object is a vector of which indices in the phylogeny tips should go first.
 This effectively ensures that the tip labels line up with the output from this function.
 It's equivalent to the following R code, where `ordered_tip_labels` is a character vector
 of the tip names in the order you always want them:
`spp_order <- match(ordered_tip_labels, phy$tip.label)`
*/
std::vector<uint> match_(std::vector<std::string> x, std::vector<std::string> table) {

    std::vector<uint> out(x.size());

    for (uint i = 0; i < out.size(); i++) {
        auto iter = std::find(table.begin(), table.end(), x[i]);
        if (iter == table.end()) stop("item in `table` not found.");
        out[i] = iter - table.begin();
    }

    return out;
}




//' Test sampling based on an evolutionary model.
//'
//' Make SURE `sampler_base_sexp` is a `ChunkMutationSampler`, not a `MutationSampler`!
//'
//' @param tip_labels Character vector of the actual phylogeny's tip labels.
//' @param ordered_tip_labels Character vector of the tip labels in the order
//'     you want them.
//'
//' @noRd
//'
//[[Rcpp::export]]
void test_phylo(SEXP& vs_sexp,
                SEXP& sampler_base_sexp,
                const uint& seq_ind,
                const std::vector<double>& branch_lens,
                const arma::Mat<uint>& edges,
                const std::vector<std::string>& tip_labels,
                const std::vector<std::string>& ordered_tip_labels,
                const arma::mat& gamma_mat) {

    XPtr<VarSet> vs_xptr(vs_sexp);
    XPtr<ChunkMutationSampler> sampler_base_xptr(sampler_base_sexp);

    std::vector<uint> spp_order = match_(ordered_tip_labels, tip_labels);

    uint n_tips = vs_xptr->size();
    if (spp_order.size() != n_tips) {
        stop("spp_order must have the same length as # variants.");
    }

    uint n_edges = edges.n_rows;
    if (branch_lens.size() != n_edges) {
        stop("branch_lens must have the same length as the # rows in edges.");
    }
    if (edges.n_cols != 2) stop("edges must have exactly two columns.");

    pcg32 eng = seeded_pcg();

    one_tree_no_recomb_<ChunkMutationSampler>(*vs_xptr, *sampler_base_xptr, seq_ind,
                                              branch_lens, edges, spp_order,
                                              gamma_mat, eng);

    return;
}
