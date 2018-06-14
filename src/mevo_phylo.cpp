
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


#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "molecular_evolution.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "mevo_phylo.h"


using namespace Rcpp;




/*
 `spp_order` should be a vector of which indices in the phylogeny tips should go first.
 This effectively ensures that the tip labels line up with the output from this function.
 In R, you could do this, where `ordered_labels` is a character vector of the tip
 names in the order you always want them:
    `spp_order <- match(ordered_labels, phy$tip.label)`
 */


void one_tree_no_recomb(VarSet& vars,
                        const uint& seq_ind,
                        const std::vector<double>& branch_lens,
                        const arma::Mat<uint>& edges,
                        const std::vector<uint>& spp_order,
                        const std::vector<std::vector<double>>& probs,
                        const std::vector<sint>& mut_lengths,
                        const std::vector<double>& pi_tcag,
                        const arma::mat& gamma_mat) {

    // arma::Mat<uint> edges = arma::conv_to<arma::Mat<uint>>::from(edges_);

    uint n_tips = vars.size();
    if (spp_order.size() != n_tips) stop("spp_order must have the same length as # variants.");

    uint n_edges = edges.n_rows;
    if (branch_lens.size() != n_edges) {
        stop("branch_lens must have the same length as the # rows in edges.");
    }
    if (edges.n_cols != 2) stop("edges must have exactly two columns.");

    uint tree_size = edges.max() + 1;

    /*
     Create tree of empty VarSequence objects
     */
    std::vector<VarSequence> var_seqs(tree_size, VarSequence(vars.reference[seq_ind]));

    /*
     Create corresponding tree of MutationSampler objects
     */
    std::vector<MutationSampler> samplers(tree_size, MutationSampler());
    for (uint i = 0; i < tree_size; i++) {
        samplers[i] = make_mutation_sampler(var_seqs[i], probs, mut_lengths,
                                            pi_tcag, gamma_mat);
    }

    // RNG
    pcg32 eng = seeded_pcg();

    /*
     `one_tree_no_recomb_` does most of the work for this function:
     */
    one_tree_no_recomb_<MutationSampler>(samplers, branch_lens, edges,
                                         n_tips, eng);

    /*
     Update final VarSequence objects in VarSet at sequence index `seq_ind`:
     */
    for (uint i = 0; i < n_tips; i++) {
        uint j = spp_order[i];
        vars[i][seq_ind].replace(var_seqs[j]);
    }

    return;

}



