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
#include <progress.hpp>  // for the progress bar


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
inline std::vector<uint32> match_(const std::vector<std::string>& ordered_tip_labels,
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
 Phylogenetic tree info for one sequence range.

 * Do all input checks outside this object's creation!
 * `edges` is converted to C++ indices inside this object's constructor.

 */
struct PhyloTree {

    std::vector<double> branch_lens;
    arma::Mat<uint32> edges;
    std::vector<std::string> tip_labels;
    uint32 start;
    sint64 end;
    std::vector<sint64> ends;  // `end` values for each tree node and tip
    uint32 n_tips;             // # tips = # variants
    uint32 tree_size;          // tree size is # tips plus # nodes
    uint32 n_edges;            // # edges = # connections between nodes/tips

    PhyloTree() {}
    PhyloTree(
        const std::vector<double>& branch_lens_,
        const arma::Mat<uint32>& edges_,
        const std::vector<std::string>& tip_labels_,
        const uint32& start_,
        const sint64& end_
    )
        : branch_lens(branch_lens_), edges(edges_), tip_labels(tip_labels_),
          start(start_), end(end_),
          ends(edges_.max(), end_),
          n_tips(tip_labels_.size()), tree_size(edges_.max()),
          n_edges(edges_.n_rows) {

        // From R to C++ indices
        edges -= 1;
    }
    PhyloTree(const PhyloTree& other)
        : branch_lens(other.branch_lens), edges(other.edges),
          tip_labels(other.tip_labels), start(other.start), end(other.end),
          ends(other.ends), n_tips(other.n_tips), tree_size(other.tree_size),
          n_edges(other.n_edges) {}

    PhyloTree& operator=(const PhyloTree& other) {
        branch_lens = other.branch_lens;
        edges = other.edges;
        tip_labels = other.tip_labels;
        start = other.start;
        end = other.end;
        ends = other.ends;
        n_tips = other.n_tips;
        tree_size = other.tree_size;
        n_edges = other.n_edges;
        return *this;
    }

};


/*
 Phylogenetic tree info for one entire sequence.

 The `trees` field contains the actual phylogeny information.
 It is a vector of `PhyloTree` objects, to allow for
 multiple phylogenies, as will be the case when recombination is included.
 If not including recombination, `trees` will simply be of length 1.

 `T` should be `MutationSampler` or `ChunkMutationSampler`.

 */


template <typename T>
class PhyloOneSeq {

public:
    std::vector<PhyloTree> trees;
    std::vector<VarSequence*> var_seq_ptrs;    // pointers to original VarSequence objects
    std::vector<VarSequence> var_seqs;  // blank VarSequence objects to evolve across tree
    std::vector<T> samplers;         // samplers to do the mutation additions across tree
    std::vector<double> seq_rates;   // sequence rates

    PhyloOneSeq() {}
    PhyloOneSeq(
        VarSet& var_set,
        const T& sampler_base,
        const uint32& seq_ind,
        const arma::mat& gamma_mat_,
        const std::vector<uint32>& n_bases_,
        const std::vector<std::vector<double>>& branch_lens_,
        const std::vector<arma::Mat<uint32>>& edges_,
        const std::vector<std::vector<std::string>>& tip_labels_
    )
        : trees(edges_.size()),
          var_seq_ptrs(var_set.size()),
          var_seqs(),
          samplers(),
          seq_rates(edges_[0].max()),
          gamma_mat(gamma_mat_),
          ordered_tip_labels(),
          recombination(branch_lens_.size() > 1)
    {

        uint32 tree_size = edges_[0].max();
        uint32 n_vars = var_set.size();
        uint32 n_trees = edges_.size();

        if (branch_lens_.size() != n_trees ||
            n_bases_.size() != n_trees ||
            tip_labels_.size() != n_trees) {
            std::string err_msg = "\nVectors for number of bases, branch lengths, ";
            err_msg += "edges, and tip labels do not all have the same length.";
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        ordered_tip_labels = std::vector<std::string>(n_vars);
        for (uint32 i = 0; i < n_vars; i++) ordered_tip_labels[i] = var_set[i].name;

        // Fill `trees`:
        uint32 start_ = 0;
        sint64 end_ = -1;
        for (uint32 i = 0; i < n_trees; i++) {
            end_ += n_bases_[i];
            trees[i] = PhyloTree(branch_lens_[i], edges_[i], tip_labels_[i],
                                 start_, end_);
            if (i < (n_trees - 1)) start_ += n_bases_[i+1];
        }

        // Filling in pointers:
        for (uint32 i = 0; i < n_vars; i++) {
            var_seq_ptrs[i] = &var_set[i][seq_ind];
        }

        // Fill in blank VarSequence objects:
        var_seqs = std::vector<VarSequence>(tree_size,
                                            VarSequence(var_set.reference[seq_ind]));

        // Fill in samplers:
        samplers = std::vector<T>(tree_size, sampler_base);
        for (uint32 i = 0; i < tree_size; i++) {
            samplers[i].fill_ptrs(var_seqs[i]);
            samplers[i].fill_gamma(gamma_mat);
        }

    }

    /*
     Similar to above, but no VarSet or sampler info yet available
     */

    PhyloOneSeq(
        const std::vector<uint32>& n_bases_,
        const std::vector<std::vector<double>>& branch_lens_,
        const std::vector<arma::Mat<uint32>>& edges_,
        const std::vector<std::vector<std::string>>& tip_labels_
    )
        : trees(edges_.size()),
          var_seq_ptrs(),
          var_seqs(),
          samplers(),
          seq_rates(edges_[0].max()),
          gamma_mat(),
          ordered_tip_labels(),
          recombination(branch_lens_.size() > 1)
    {

        uint32 n_trees = edges_.size();

        if (n_bases_.size() != branch_lens_.size() ||
            n_bases_.size() != edges_.size() ||
            n_bases_.size() != tip_labels_.size()) {
            std::string err_msg = "\nVectors for number of bases, branch lengths, edges";
            err_msg += ", and tip labels do not all have the same length.";
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        // Fill `trees`:
        uint32 start_ = 0;
        sint64 end_ = -1;
        for (uint32 i = 0; i < n_trees; i++) {
            end_ += n_bases_[i];
            trees[i] = PhyloTree(branch_lens_[i], edges_[i], tip_labels_[i],
                                 start_, end_);
            if (i < (n_trees - 1)) start_ += n_bases_[i+1];
        }

    }

    /*
     Set sampler and variant info:
     */
    void set_samp_var_info(VarSet& var_set,
                           const T& sampler_base,
                           const uint32& seq_ind,
                           const arma::mat& gamma_mat_) {

        uint32 tree_size = trees[0].tree_size;
        uint32 n_vars = var_set.size();

        var_seq_ptrs.resize(var_set.size());
        gamma_mat = gamma_mat_;

        ordered_tip_labels = std::vector<std::string>(n_vars);
        for (uint32 i = 0; i < n_vars; i++) ordered_tip_labels[i] = var_set[i].name;

        // Filling in pointers:
        for (uint32 i = 0; i < n_vars; i++) {
            var_seq_ptrs[i] = &var_set[i][seq_ind];
        }

        // Fill in blank VarSequence objects:
        var_seqs = std::vector<VarSequence>(tree_size,
                                            VarSequence(var_set.reference[seq_ind]));

        // Fill in samplers:
        samplers = std::vector<T>(tree_size, sampler_base);
        for (uint32 i = 0; i < tree_size; i++) {
            samplers[i].fill_ptrs(var_seqs[i]);
            samplers[i].fill_gamma(gamma_mat);
        }

        return;

    }



    /*
     Evolve all trees.
     */
    int evolve(pcg32& eng, Progress& prog_bar) {

        for (PhyloTree& tree : trees) {
            int code = one_tree(tree, eng, prog_bar);
            if (code == -1) return code;
        }
        return 0;
    }



private:
    arma::mat gamma_mat;
    std::vector<std::string> ordered_tip_labels;
    bool recombination;


    /*
     Evolve one tree.
     */
    int one_tree(PhyloTree& tree, pcg32& eng, Progress& prog_bar);


    /*
     Reset for a new tree:
     */
    void reset(const PhyloTree& tree) {

        const uint32& tree_size(tree.tree_size);
        const uint32& start(tree.start);
        const uint32& end(tree.end);

        if (tree_size == 0) {
            throw(Rcpp::exception("\ntree size of zero is non-sensical.", false));
        }
        // Resize blank VarSequence objects if necessary:
        if (tree_size != var_seqs.size()) {
            VarSequence var_seq_(var_seq_ptrs[0]->ref_seq);
            var_seqs.resize(tree_size, var_seq_);
        }
        // Empty mutations from tree of VarSequence objects:
        for (VarSequence& var_seq : var_seqs) var_seq.clear();

        // Re-size samplers if necessary
        if (tree_size != samplers.size()) {
            samplers.resize(tree_size, samplers[0]);
        }
        // Fill in sampler pointers and original gamma matrix:
        for (uint32 i = 0; i < tree_size; i++) {
            samplers[i].fill_ptrs(var_seqs[i]);
            samplers[i].fill_gamma(gamma_mat);
        }
        /*
         Set up vector of overall sequence rates.
         They should all be the same as the first one to start out.
         */
        double rate_;
        if (!recombination) {
            rate_ = samplers[0].total_rate();
        } else rate_ = samplers[0].total_rate(start, static_cast<uint32>(end), true);
        seq_rates = std::vector<double>(tree_size, rate_);

        return;
    }


    /*
     Update for a new edge:
     */
    void update(std::exponential_distribution<double>& distr,
                const uint32& b1,
                const uint32& b2) {

        /*
         Replace existing mutation information in VarSequence at `b1` with info in the
         one at `b2`
         */
        samplers[b2].var_seq->replace(*samplers[b1].var_seq);
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
     Clear info from VarSequence object at `b1` if it's no longer needed, to free up
     some memory.
     (If it's the last branch length, `b1` will always be a node and thus no longer
     needed.)
     */
    void clear_branches(const uint32& b1,
                        const uint32& i,
                        const PhyloTree& tree) {
        bool clear_b1;
        if (i < (tree.edges.n_rows - 1)) {
            // Is it absent from any remaining items in the first column?
            clear_b1 = ! arma::any(tree.edges(arma::span(i+1, tree.edges.n_rows - 1),
                                              0) == b1);
        } else clear_b1 = true;
        if (clear_b1) samplers[b1].var_seq->clear();
        return;
    }

    /*
     Update final `VarSequence` objects in the same order as `ordered_tip_labels`:
     */
    void fill(const PhyloTree& tree) {

        std::vector<uint32> spp_order = match_(ordered_tip_labels,
                                               tree.tip_labels);

        if (recombination) {
            for (uint32 i = 0; i < tree.n_tips; i++) {
                uint32 j = spp_order[i];
                (*var_seq_ptrs[i]) += var_seqs[j];
            }
        } else {
            for (uint32 i = 0; i < tree.n_tips; i++) {
                uint32 j = spp_order[i];
                (*var_seq_ptrs[i]).replace(var_seqs[j]);
            }
        }
        return;
    }

};



/*
 Process one phylogenetic tree for a single sequence with no recombination.
 This template does most of the work for the chunked and non-chunked versions in
 the cpp file.
 `T` should be `MutationSampler` or `ChunkMutationSampler`.

 Note that this function should be changed if any of these VarSequences differ from
 each other (within the range specified if recombination = true).
 They can already have mutations, but to start out, they must all be the same.
 */
template<typename T>
int PhyloOneSeq<T>::one_tree(PhyloTree& tree,
                             pcg32& eng,
                             Progress& prog_bar) {

    // Reset tree of samplers and VarSequence objects representing nodes and tips:
    reset(tree);

    /*
     Check for a user interrupt. Using a Progress object allows the user to interrupt
     the process during multithreaded operations.
     If recombination == true, I'm only doing this here, not for each edge bc that
     would likely cause too many checks, which would slow things down.
     */
    if (prog_bar.is_aborted()) return -1;

    // Exponential distribution to do the time-jumps along the branch lengths:
    std::exponential_distribution<double> distr(1.0);

    /*
     Now iterate through the phylogeny:
     */
    for (uint32 i = 0; i < tree.n_edges; i++) {

        // If not simulating recombination, checking for abort every edge:
        if (!recombination) {
            if (prog_bar.is_aborted()) return -1;
        }

        // Indices for nodes/tips that the branch length in `branch_lens` refers to
        uint32 b1 = tree.edges(i,0);
        uint32 b2 = tree.edges(i,1);

        /*
         Update `samplers`, `seq_rates`, and `distr` for this edge:
         */
        update(distr, b1, b2);

        /*
         Now do exponential jumps and mutate until you exceed the branch length.
         */
        double& rate(seq_rates[b2]);
        double amt_time = tree.branch_lens[i];
        double time_jumped = distr(eng);
        double rate_change = 0;
        if (recombination) {
            sint64& end_(tree.ends[b2]);
            end_ = tree.ends[b1];
            const sint64 start_ = static_cast<sint64>(tree.start);
            while (time_jumped <= amt_time && end_ >= start_) {
                /*
                 Add mutation here, outputting how much the overall sequence rate should
                 change:
                 (`end_` is automatically adjusted for indels)
                 */
                rate_change = samplers[b2].mutate(eng, tree.start, end_);
                /*
                 Adjust the overall sequence rate, then update the exponential
                 distribution:
                 */
                rate += rate_change;
                distr.param(std::exponential_distribution<double>::param_type(rate));
                // Jump again:
                time_jumped += distr(eng);
            }
        } else {
            // Same thing but without recombination
            while (time_jumped <= amt_time && var_seqs[b2].size() > 0) {
                rate_change = samplers[b2].mutate(eng);
                rate += rate_change;
                distr.param(std::exponential_distribution<double>::param_type(rate));
                time_jumped += distr(eng);
            }
        }

        /*
         To free up some memory, clear info from VarSequence object at `b1` if it's no
         longer needed.
         */
        clear_branches(b1, i, tree);

    }

    /*
     Update final `VarSequence` objects:
     */
    fill(tree);

    // Update progress bar:
    if (recombination) {
        prog_bar.increment(tree.end - tree.start + 1);
    } else prog_bar.increment(var_seq_ptrs[0]->ref_seq.size());

    return 0;

}







/*
 `T` should be `MutationSampler` or `ChunkMutationSampler`.
 */
template <typename T>
void evolve_seqs_(
        SEXP& var_set_ptr,
        SEXP& sampler_base_ptr,
        SEXP& phylo_info_ptr,
        const std::vector<uint32>& seq_inds,
        const std::vector<arma::mat>& gamma_mats,
        const bool& show_progress) {

    XPtr<VarSet> var_set(var_set_ptr);
    XPtr<T> sampler_base(sampler_base_ptr);
    XPtr<std::vector<PhyloOneSeq<T>>> phylo_info(phylo_info_ptr);

    uint32 n_seqs = var_set->reference.size();
    uint64 total_seq = var_set->reference.total_size;

    Progress prog_bar(total_seq, show_progress);

    if (n_seqs != gamma_mats.size()) {
        std::string err_msg = "\ngamma_mats must be of same length as # sequences in ";
        err_msg += "reference";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }
    if (n_seqs != seq_inds.size()) {
        std::string err_msg = "\nseq_inds must be of same length as # sequences in ";
        err_msg += "reference";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }
    if (n_seqs != phylo_info->size()) {
        std::string err_msg = "\nphylo_info must be of same length as # sequences in ";
        err_msg += "reference";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }


    pcg32 eng = seeded_pcg();

    for (uint32 i = 0; i < n_seqs; i++) {

        PhyloOneSeq<T>& seq_phylo((*phylo_info)[i]);

        const arma::mat& gamma_mat(gamma_mats[i]);
        const uint32& seq_ind(seq_inds[i]);

        if (gamma_mat(gamma_mat.n_rows-1,0) != (*var_set)[0][seq_ind].size()) {
            std::string err_msg = "\nGamma matrices must have max values equal to ";
            err_msg += "the respective sequence's length.\n";
            err_msg += "This error occurred on Gamma matrix number ";
            err_msg += std::to_string(i+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        // Set values for variant info and sampler:
        seq_phylo.set_samp_var_info(*var_set, *sampler_base, seq_ind, gamma_mat);

        // Evolve the sequence using the seq_phylo object:
        int code = seq_phylo.evolve(eng, prog_bar);

        // If the code is -1, the user interrupted the process.
        // Make sure this check happens outside of multithreaded code.
        if (code == -1) {
            std::string warn_msg = "\nThe user interrupted phylogenetic evolution. ";
            warn_msg += "Note that changes occur in place, so your variants have ";
            warn_msg += "already been partially added.";
            Rcpp::warning(warn_msg.c_str());
            return;
        }
    }

    return;
}






#endif
