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


using namespace Rcpp;




//' Process one phylogenetic tree for a single sequence with no recombination.
//'
//' @noRd
//'
void one_tree_no_recomb(const VarSet& vars,
                        const uint& seq_ind,
                        double& seq_rate,
                        const std::vector<std::string>& labels,
                        const std::vector<double>& branch_lens,
                        const arma::mat& edges) {

    uint n_tips = vars.size();
    if (labels.size() != n_tips) stop("labels must have the same length as # variants.");
    if (edges.n_rows != n_tips) stop("edges must have the same nrows as # variants.");
    if (branch_lens.size() != n_tips) {
        stop("branch_lens must have the same length as # variants.");
    }

    /*
     Create tree of VarSequence objects
     */
    ;
    /*
     Create corresponding tree of MutationSampler or ChunkMutationSampler objects
     */
    ;


    /*
     Now iterate through phylogeny:
     */
    for (uint i = 0; i < n_tips; i++) {
        uint b1 = static_cast<uint>(edges(i,1)) - 1;
        uint b2 = static_cast<uint>(edges(i,2)) - 1;
        uint muts = 0;
        double time_left = branch_lens[i];
        double t_jump = R::rexp(seq_rate);
        while (t_jump <= time_left) {
            /*
             Add mutation here
             */
            muts++;
            /*
             Adjust `seq_rate`
             */
            t_jump += R::rexp(seq_rate);
        }
        bool clear_b1;
        if (i < (n_tips - 1)) {
            clear_b1 = arma::any(edges(arma::span(i+1, edges.n_rows - 1), 0) == b1);
        } else clear_b1 = true;
        if (clear_b1) {
            /*
             Erase VarSequence object at `b1`
             */
            ;
        }
        if (b2 < n_tips) {
            /*
             Add final VarSequence object to VarSet at index `b2`
             */
            ;
        }
    }


    /*
# Will be creating tree of VarGenome objects, NOT VarSet objects, so don't worry about
# copying RefGenome objects. This won't be a problem.

    for (i in 1:nrow(edges)) {
    b1 <- edges[i,1]
    b2 <- edges[i,2]
    muts <- 0
    time_left <- branch_lens[i]
    t_jump <- rexp(1, rate = seq_rate)
    while (t_jump < time_left) {
    muts <- muts+1
    t_jump = t_jump + rexp(1, rate = seq_rate)
    }
    cat(sprintf("%i mutations from %i to %i\n", muts, b1, b2))
    if (i < nrow(edges)) {
    if (!b1 %in% edges[(i+1):(nrow(edges)),1]) cat(sprintf("remove %i at %i\n", b1, i))
    } else cat(sprintf("remove %i at %i\n", b1, i))
    if (b2 <= n_tips) cat(sprintf("%i added at %i\n", b2, i))
    }

     */

}


#endif
