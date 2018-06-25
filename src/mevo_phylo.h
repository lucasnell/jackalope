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
#include "mevo.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types


using namespace Rcpp;





/*
This function produces a vector of which indices in the phylogeny tips should go first.
This effectively ensures that the tip labels line up with the output from this function.
It's equivalent to the following R code, where `ordered_tip_labels` is a character vector
of the tip names in the order you always want them:
`spp_order <- match(ordered_tip_labels, phy$tip.label)`
*/
std::vector<uint32> match_(const std::vector<std::string>& ordered_tip_labels,
                           const std::vector<std::string>& tip_labels) {

    std::vector<uint32> spp_order(ordered_tip_labels.size());

    for (uint32 i = 0; i < spp_order.size(); i++) {
        auto iter = std::find(tip_labels.begin(), tip_labels.end(),
                              ordered_tip_labels[i]);
        if (iter == tip_labels.end()) stop("item in `tip_labels` not found.");
        spp_order[i] = iter - tip_labels.begin();
    }

    return spp_order;
}




/*
 Update `samplers`, `seq_rates`, and `distr` for a given edge.

 <Recombination>
 This works the same with or without recombination.
*/
template <typename T>
void update_samplers_rates_distr_(std::vector<T>& samplers,
                                  std::vector<double>& seq_rates,
                                  std::exponential_distribution<double>& distr,
                                  const uint32& b1, const uint32& b2) {

    /*
     Replace existing mutation information in VarSequence at `b1` with info in the
     one at `b2`
     */
    samplers[b2].vs->replace(*samplers[b1].vs);
    /*
     Do the same for the SeqGammas in the sampler:
     */
    samplers[b2].location.mr().gammas = samplers[b1].location.mr().gammas;
    /*
     Update overall sequence rate:
     */
    seq_rates[b2] = seq_rates[b1];

    // Set exponential distribution to use this sequence's rate:
    distr.param(std::exponential_distribution<double>::param_type(seq_rates[b2]));

    return;
}


/*
 Create tree of empty VarSequence objects, corresponding tree of
 [Chunk]MutationSampler objects, and fill in sequence rates.

 <Recombination>
 To do this for a whole sequence (i.e., without recombination), just
 use default values for `whole_sequence`, `start`, and `end`.
 If doing this for a region (i.e., with recombination), make sure `start` and `end`
 are specified and that `whole_sequence` is false!

 */
template <typename T>
void fill_var_samp_rate_(std::vector<VarSequence>& var_seqs,
                         std::vector<T>& samplers,
                         std::vector<double>& seq_rates,
                         const uint32& tree_size, const VarSet& vars,
                         const uint32& seq_ind, const T& sampler_base,
                         const arma::mat& gamma_mat,
                         const bool& whole_sequence = true,
                         const uint32& start = 0, const uint32& end = 0) {

    var_seqs = std::vector<VarSequence>(tree_size, VarSequence(vars.reference[seq_ind]));

    /*
     For tree of [Chunk]MutationSampler objects,
     first fill with the base sampler, then fill in pointers to the corresponding
     VarSequence object, then fill in the gamma matrix.
     */
    samplers = std::vector<T>(tree_size, sampler_base);
    for (uint32 i = 0; i < tree_size; i++) {
        samplers[i].fill_ptrs(var_seqs[i]);
        samplers[i].fill_gamma(gamma_mat);
    }
    /*
     Set up vector of overall sequence rates.
     They should all be the same as the first one to start out.
     */
    double rate_;
    if (whole_sequence) {
        rate_ = samplers[0].total_rate();
    } else rate_ = samplers[0].total_rate(start, end, true);
    seq_rates = std::vector<double>(tree_size, rate_);

    return;
}

/*
 Clear info from VarSequence object at `b1` if it's no longer needed to
 free up some memory.
 (If it's the last branch length, `b1` will always be a node and thus no longer
 needed.)
 */
template <typename T>
void clear_branches_(std::vector<T>& samplers,
                     const uint32& b1,
                     const uint32& i,
                     const arma::Mat<uint32>& edges) {
    bool clear_b1;
    if (i < (edges.n_rows - 1)) {
        // Is it absent from any remaining items in the first column?
        clear_b1 = ! arma::any(edges(arma::span(i+1, edges.n_rows - 1), 0) == b1);
    } else clear_b1 = true;
    if (clear_b1) samplers[b1].vs->clear();
    return;
}







#endif
