#ifndef __GEMINO_MEVO_PHYLO_H
#define __GEMINO_MEVO_PHYLO_H


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
#include <random>  // exponential_distribution


#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "molecular_evolution.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types


using namespace Rcpp;



/*
 Process one phylogenetic tree for a single sequence with no recombination.
 This template does most of the work for the chunked and non-chunked versions in
 the cpp file.
 `T` should be `MutationSampler` or `ChunkMutationSampler`.

 Note that this function should be changed if any of these VarSequences differ from
 each other.
 They can already have mutations, but to start out, they must all be the same.

 Also, this function does not do any reordering of branches in relation to species names.
 In other words, the vector of VarSequence objects is indexed by tip numbers.
 */
template <typename T>
inline void one_tree_no_recomb_(VarSet& vars,
                                const T& sampler_base,
                                const uint& seq_ind,
                                const std::vector<double>& branch_lens,
                                const arma::Mat<uint>& edges,
                                const std::vector<uint>& spp_order,
                                const arma::mat& gamma_mat,
                                pcg32& eng) {

    // # tips = # variants
    uint n_tips = vars.size();
    // tree size is # tips plus # nodes
    uint tree_size = edges.max() + 1;
    // Number of edges = the number of connections between nodes/tips
    uint n_edges = edges.n_rows;

    /*
     Create tree of the same VarSequence objects
     */
    std::vector<VarSequence*> var_seqs(tree_size);
    for (uint i = 0; i < tree_size; i++) var_seqs[i] = new VarSequence(vars.reference[seq_ind]);

    /*
     Create corresponding tree of MutationSampler objects
     First fill with the base sampler, then fill in pointers to the corresponding
     VarSequence object, then fill in the gamma matrix.
     */
    std::vector<T> samplers(tree_size, sampler_base);
    for (uint i = 0; i < tree_size; i++) {
        samplers[i].fill_ptrs(*var_seqs[i]);
        samplers[i].fill_gamma(gamma_mat);
    }


    /*
     Exponential distribution to do the time-jumps along the branch lengths:
     */
    std::exponential_distribution<double> distr;

    /*
     Set up vector of overall sequence rates.
     They should all be the same as the first one to start out.
     */
    std::vector<double> seq_rates(tree_size, samplers[0].total_rate());

    /*
     Now iterate through the phylogeny:
     */
    for (uint i = 0; i < n_edges; i++) {
        // Indices for nodes/tips that the branch length in `branch_lens` refers to
        uint b1 = edges(i,0);
        uint b2 = edges(i,1);
        /*
         Replace existing mutation information in VarSequence at `b1` with info in the
         one at `b2`
         */
        samplers[b2].vs->replace(*samplers[b1].vs);
        /*
         Update overall sequence rate:
         */
        double& rate(seq_rates[b2]);
        rate = seq_rates[b1];

        // Set exponential distribution to use this sequence's rate:
        distr.param(std::exponential_distribution<double>::param_type(rate));

        /*
         Now do exponential jumps and mutate until you exceed the branch length.
         */
        double amt_time = branch_lens[i];
        double time_jumped = distr(eng);
        while (time_jumped <= amt_time) {
            /*
             Add mutation here, outputting how much the overall sequence rate should
             change:
             */
            double rate_change = samplers[b2].mutate_rate_change(eng);
            /*
             Adjust the overall sequence rate, then update the exponential distribution:
             */
            rate += rate_change;
            distr.param(std::exponential_distribution<double>::param_type(rate));
            /*
             Jump again:
             */
            time_jumped += distr(eng);
        }
        /*
         Clear info from VarSequence object at `b1` if it's no longer needed to
         free up some memory.
         (If it's the last branch length, `b1` will always be a node and thus no longer
         needed.)
         */
        bool clear_b1;
        if (i < (n_tips - 1)) {
            // Is it absent from any remaining items in the first column?
            clear_b1 = ! arma::any(edges(arma::span(i+1, edges.n_rows - 1), 0) == b1);
        } else clear_b1 = true;
        if (clear_b1) samplers[b1].vs->clear();
    }

    /*
     Update final `VarSequence` objects using the index vector `spp_order`:
     */
    for (uint i = 0; i < n_tips; i++) {
        uint j = spp_order[i];
        vars[i][seq_ind].replace(*var_seqs[j]);
    }

    return;

}









#endif
