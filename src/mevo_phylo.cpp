// Below is added to entirely prevent this from compiling for now
#define __GEMINO_MEVO_PHYLO_CPP


#ifndef __GEMINO_MEVO_PHYLO_CPP
#define __GEMINO_MEVO_PHYLO_CPP




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
 */


void one_tree_no_recomb(VarSet& vars,
                        const uint& seq_ind,
                        const std::vector<double>& branch_lens,
                        const arma::mat& edges_,
                        const std::vector<uint>& spp_order,
                        const std::vector<std::vector<double>>& probs,
                        const std::vector<sint>& mut_lengths,
                        const std::vector<double>& pi_tcag,
                        const arma::mat& gamma_mat) {

    arma::Mat<uint> edges = arma::conv_to<arma::Mat<uint>>::from(edges_);

    uint n_tips = vars.size();
    if (spp_order.size() != n_tips) stop("spp_order must have the same length as # variants.");

    uint n_edges = edges.n_rows;
    if (branch_lens.size() != n_edges) {
        stop("branch_lens must have the same length as the # rows in edges.");
    }
    if (edges.n_cols != 2) stop("edges must have exactly two columns.");

    uint tree_size = edges.max();


    /*
    Create tree of empty VarSequence objects
    */
    std::vector<VarSequence> var_seqs(tree_size, VarSequence(vars.reference[seq_ind]));

    /*
     Create corresponding tree of MutationSampler or ChunkMutationSampler objects
     */
    std::vector<MutationSampler> samplers;
    for (uint i = 0; i < tree_size; i++) {
        MutationSampler ms = make_mutation_sampler(var_seqs[i], probs, mut_lengths,
                                                   pi_tcag, gamma_mat);
        samplers.push_back(ms);
    }

    // RNG
    pcg32 eng = seeded_pcg();

    /*
     Does most of the work for this function:
     */
    one_tree_no_recomb_<MutationSampler>(var_seqs, samplers, branch_lens, edges, eng);

    /*
     Update final VarSequence objects in VarSet at sequence index `seq_ind`:
     */
    for (uint i = 0; i < n_tips; i++) {
        uint j = spp_order[i];
        vars[i][seq_ind].mutations = var_seqs[j].mutations;
        vars[i][seq_ind].seq_size = var_seqs[j].seq_size;
    }

    return;

}


#endif
