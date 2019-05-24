#ifndef __JACKAL_PHYLOGENOMICS_H
#define __JACKAL_PHYLOGENOMICS_H


/*
 ********************************************************

 Methods for evolving sequences along phylogenies / gene trees

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <random>  // exponential_distribution
#include <progress.hpp>  // for the progress bar
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // integer types
#include "seq_classes_var.h"  // Var* classes
#include "mutator.h"  // samplers
#include "alias_sampler.h" // alias sampling
#include "pcg.h" // pcg sampler types


using namespace Rcpp;





/*
 This function produces a vector of which indices in the phylogeny tips should go first.
 This effectively ensures that the tip labels line up with the output from this function.
 It's equivalent to the following R code, where `ordered_tip_labels` is a character vector
 of the tip names in the order you always want them:
 `spp_order <- match(ordered_tip_labels, phy$tip.label)`
 */
inline std::vector<uint64> match_(const std::vector<std::string>& ordered_tip_labels,
                                  const std::vector<std::string>& tip_labels) {

    std::vector<uint64> spp_order(ordered_tip_labels.size());

    for (uint64 i = 0; i < spp_order.size(); i++) {
        auto iter = std::find(tip_labels.begin(), tip_labels.end(),
                              ordered_tip_labels[i]);
        if (iter == tip_labels.end()) {
            std::string err_msg = ordered_tip_labels[i];
            err_msg += " not found in `tip_labels`.";
            stop(err_msg.c_str());
        }
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
    arma::Mat<uint64> edges;
    std::vector<std::string> tip_labels;
    uint64 start;
    sint64 end;
    std::vector<sint64> ends;  // `end` values for each tree node and tip
    uint64 n_tips;             // # tips = # variants
    uint64 tree_size;          // tree size is # tips plus # nodes
    uint64 n_edges;            // # edges = # connections between nodes/tips

    PhyloTree() {}
    PhyloTree(
        const std::vector<double>& branch_lens_,
        const arma::Mat<uint64>& edges_,
        const std::vector<std::string>& tip_labels_,
        const uint64& start_,
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

 */


class PhyloOneSeq {

public:
    std::vector<PhyloTree> trees;
    std::vector<VarSequence*> var_seq_ptrs;    // pointers to original VarSequence objects
    std::vector<VarSequence> var_seqs;  // blank VarSequence objects to evolve across tree
    std::vector<MutationSampler> samplers; // to do the mutation additions across tree
    std::vector<double> seq_rates;   // sequence rates
    uint64 n_tips;                  // number of tips (i.e., variants)

    PhyloOneSeq() {}
    PhyloOneSeq(
        VarSet& var_set,
        const MutationSampler& sampler_base,
        const uint64& seq_ind,
        const arma::mat& gamma_mat_,
        const std::vector<uint64>& n_bases_,
        const std::vector<std::vector<double>>& branch_lens_,
        const std::vector<arma::Mat<uint64>>& edges_,
        const std::vector<std::vector<std::string>>& tip_labels_
    )
        : trees(edges_.size()),
          var_seq_ptrs(var_set.size()),
          var_seqs(),
          samplers(),
          seq_rates(edges_[0].max()),
          n_tips(tip_labels_[0].size()),
          gamma_mat(gamma_mat_),
          ordered_tip_labels(),
          recombination(branch_lens_.size() > 1)
    {

        uint64 tree_size = edges_[0].max();
        uint64 n_vars = var_set.size();
        uint64 n_trees = edges_.size();

        if (branch_lens_.size() != n_trees ||
            n_bases_.size() != n_trees ||
            tip_labels_.size() != n_trees) {
            std::string err_msg = "\nVectors for number of bases, branch lengths, ";
            err_msg += "edges, and tip labels do not all have the same length.";
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        ordered_tip_labels = std::vector<std::string>(n_vars);
        for (uint64 i = 0; i < n_vars; i++) ordered_tip_labels[i] = var_set[i].name;

        // Fill `trees`:
        uint64 start_ = 0;
        sint64 end_ = -1;
        for (uint64 i = 0; i < n_trees; i++) {
            end_ += n_bases_[i];
            if (i > 0) start_ += n_bases_[i-1];
            trees[i] = PhyloTree(branch_lens_[i], edges_[i], tip_labels_[i],
                                 start_, end_);
        }

        // Filling in pointers:
        for (uint64 i = 0; i < n_vars; i++) {
            var_seq_ptrs[i] = &var_set[i][seq_ind];
        }

        // Fill in blank VarSequence objects:
        var_seqs = std::vector<VarSequence>(tree_size,
                                            VarSequence((*var_set.reference)[seq_ind]));

        // Fill in samplers:
        samplers = std::vector<MutationSampler>(tree_size, sampler_base);
        for (uint64 i = 0; i < tree_size; i++) {
            samplers[i].new_seq(var_seqs[i], gamma_mat);
        }

    }

    /*
     Similar to above, but no VarSet or sampler info yet available
     */

    PhyloOneSeq(
        const std::vector<uint64>& n_bases_,
        const std::vector<std::vector<double>>& branch_lens_,
        const std::vector<arma::Mat<uint64>>& edges_,
        const std::vector<std::vector<std::string>>& tip_labels_
    )
        : trees(edges_.size()),
          var_seq_ptrs(),
          var_seqs(),
          samplers(),
          seq_rates(edges_[0].max()),
          n_tips(tip_labels_[0].size()),
          gamma_mat(),
          ordered_tip_labels(),
          recombination(branch_lens_.size() > 1)
    {

        uint64 n_trees = edges_.size();

        if (n_bases_.size() != branch_lens_.size() ||
            n_bases_.size() != edges_.size() ||
            n_bases_.size() != tip_labels_.size()) {
            std::string err_msg = "\nVectors for number of bases, branch lengths, edges";
            err_msg += ", and tip labels do not all have the same length.";
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        // Fill `trees`:
        uint64 start_ = 0;
        sint64 end_ = -1;
        for (uint64 i = 0; i < n_trees; i++) {
            end_ += n_bases_[i];
            if (i > 0) start_ += n_bases_[i-1];
            trees[i] = PhyloTree(branch_lens_[i], edges_[i], tip_labels_[i],
                                 start_, end_);
        }

    }

    /*
     Set sampler and variant info:
     */
    void set_samp_var_info(VarSet& var_set,
                           const MutationSampler& sampler_base,
                           const uint64& seq_ind,
                           const arma::mat& gamma_mat_) {

        uint64 tree_size = trees[0].tree_size;
        uint64 n_vars = var_set.size();

        gamma_mat = gamma_mat_;

        ordered_tip_labels.resize(n_vars);
        for (uint64 i = 0; i < n_vars; i++) ordered_tip_labels[i] = var_set[i].name;

        // Filling in pointers:
        var_seq_ptrs.resize(n_vars);
        for (uint64 i = 0; i < n_vars; i++) {
            var_seq_ptrs[i] = &var_set[i][seq_ind];
        }

        // Fill in blank VarSequence objects:
        var_seqs = std::vector<VarSequence>(tree_size,
                                            VarSequence((*var_set.reference)[seq_ind]));

        // Fill in samplers:
        samplers = std::vector<MutationSampler>(tree_size, sampler_base);
        for (uint64 i = 0; i < tree_size; i++) {
            samplers[i].new_seq(var_seqs[i], gamma_mat);
        }

        return;

    }


    /*
     Evolve all trees.
     */
    int evolve(pcg64& eng, Progress& prog_bar) {

        for (PhyloTree& tree : trees) {
            int code = one_tree(tree, eng, prog_bar);
            if (code == -1) return code;
        }
        return 0;
    }



    /*
     Fill a PhyloOneSeq object from an input list
    */
    void fill_from_list(const List& genome_phylo_info, const uint64& i) {

        std::string err_msg;

        const List& seq_phylo_info(genome_phylo_info[i]);

        uint64 n_trees = seq_phylo_info.size();
        if (n_trees == 0) {
            err_msg = "\nNo trees supplied on sequence " + std::to_string(i+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        std::vector<uint64> n_bases_(n_trees);
        std::vector<std::vector<double>> branch_lens_(n_trees);
        std::vector<arma::Mat<uint64>> edges_(n_trees);
        std::vector<std::vector<std::string>> tip_labels_(n_trees);

        for (uint64 j = 0; j < n_trees; j++) {

            const List& phylo_info(seq_phylo_info[j]);
            std::vector<double> branch_lens = as<std::vector<double>>(
                phylo_info["branch_lens"]);
            arma::Mat<uint64> edges = as<arma::Mat<uint64>>(phylo_info["edges"]);
            std::vector<std::string> tip_labels = as<std::vector<std::string>>(
                phylo_info["labels"]);
            if (branch_lens.size() != edges.n_rows) {
                err_msg = "\nBranch lengths and edges don't have the same ";
                err_msg += "size on sequence " + std::to_string(i+1) + " and tree ";
                err_msg += std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }
            if (branch_lens.size() == 0) {
                err_msg = "\nEmpty tree on sequence " + std::to_string(i+1);
                err_msg += " and tree " + std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }
            uint64 start = as<uint64>(phylo_info["start"]);
            sint64 end = as<sint64>(phylo_info["end"]);
            if (end < static_cast<sint64>(start)) {
                err_msg = "\nEnd position < start position on sequence ";
                err_msg += std::to_string(i+1) + " and tree " + std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }

            n_bases_[j] = end - start + 1;
            branch_lens_[j] = branch_lens;
            edges_[j] = edges;
            tip_labels_[j] = tip_labels;
        }


        *this = PhyloOneSeq(n_bases_, branch_lens_, edges_, tip_labels_);

        return;

    }




private:
    arma::mat gamma_mat;
    std::vector<std::string> ordered_tip_labels;
    bool recombination;


    /*
     Evolve one tree.
     */
    int one_tree(PhyloTree& tree, pcg64& eng, Progress& prog_bar);


    /*
     Reset for a new tree:
     */
    void reset(const PhyloTree& tree) {

        const uint64& tree_size(tree.tree_size);
        const uint64& start(tree.start);
        const uint64& end(tree.end);

        if (tree_size == 0) {
            throw(Rcpp::exception("\ntree size of zero is non-sensical.", false));
        }
        // Resize blank VarSequence objects if necessary:
        if (tree_size != var_seqs.size()) {
            VarSequence var_seq_(*(var_seq_ptrs[0]->ref_seq));
            var_seqs.resize(tree_size, var_seq_);
        }
        // Empty mutations from tree of VarSequence objects:
        for (VarSequence& var_seq : var_seqs) var_seq.clear();

        // Re-size samplers if necessary
        if (tree_size != samplers.size()) {
            samplers.resize(tree_size, samplers[0]);
        }
        // Fill in sampler pointers and original gamma matrix:
        for (uint64 i = 0; i < tree_size; i++) {
            samplers[i].new_seq(var_seqs[i], gamma_mat);
        }
        /*
         Set up vector of overall sequence rates.
         They should all be the same as the first one to start out.
         */
        double rate_;
        if (!recombination) {
            rate_ = samplers[0].location.total_rate();
        } else {
            // (This will be checked every call to `sample`, so no need to do it for
            // all samplers)
            samplers[0].location.new_bounds(start, end);
            rate_ = samplers[0].location.bounds.end_rate -
                samplers[0].location.bounds.start_rate;
        }
        seq_rates = std::vector<double>(tree_size, rate_);

        return;
    }


    /*
     Update for a new edge:
     */
    void update(std::exponential_distribution<double>& distr,
                const uint64& b1,
                const uint64& b2) {

        /*
         Replace existing mutation information in VarSequence at `b1` with info in the
         one at `b2`
         */
        samplers[b2].var_seq->replace(*samplers[b1].var_seq);
        /*
         Do the same for the SeqGammas in the sampler:
         */
        samplers[b2].location.regions = samplers[b1].location.regions;

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
    void clear_branches(const uint64& b1,
                        const uint64& i,
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
    void update_var_seq(const PhyloTree& tree);

};







/*
 Phylogenetic tree info for all sequences in a genome.

 */
class PhyloInfo {
public:

    std::vector<PhyloOneSeq> phylo_one_seqs;

    PhyloInfo(const List& genome_phylo_info) {

        uint64 n_seqs = genome_phylo_info.size();

        if (n_seqs == 0) {
            throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                                  false));
        }

        phylo_one_seqs = std::vector<PhyloOneSeq>(n_seqs);

        for (uint64 i = 0; i < n_seqs; i++) {
            phylo_one_seqs[i].fill_from_list(genome_phylo_info, i);
        }
    }

    XPtr<VarSet> evolve_seqs(
            SEXP& ref_genome_ptr,
            SEXP& sampler_base_ptr,
            const std::vector<arma::mat>& gamma_mats,
            uint64 n_threads,
            const bool& show_progress);



};








#endif
