#ifndef __JACKALOPE_PHYLOGENOMICS_H
#define __JACKALOPE_PHYLOGENOMICS_H


/*
 ********************************************************

 Methods for evolving chromosomes along phylogenies / gene trees

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
#include "var_classes.h"  // Var* classes
#include "mutator.h"  // TreeMutator
#include "alias_sampler.h" // alias sampling
#include "pcg.h" // pcg sampler types


using namespace Rcpp;






/*
 Phylogenetic tree info for one chromosome range.

 * Do all input checks outside this object's creation!
 * `edges` is converted to C++ indices inside this object's constructor.

 */
struct PhyloTree {

    std::vector<double> branch_lens;
    arma::Mat<uint64> edges;
    std::vector<std::string> tip_labels;
    uint64 start;
    uint64 end;
    std::vector<uint64> ends;  // (non-inclusive) `end` values for each tree node and tip
    uint64 n_tips;             // # tips = # variants
    uint64 tree_size;          // tree size is # tips plus # nodes
    uint64 n_edges;            // # edges = # connections between nodes/tips

    PhyloTree() {}
    PhyloTree(
        const std::vector<double>& branch_lens_,
        const arma::Mat<uint64>& edges_,
        const std::vector<std::string>& tip_labels_,
        const uint64& start_,
        const uint64& end_
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
 Phylogenetic tree info for one entire chromosome.

 The `trees` field contains the actual phylogeny information.
 It is a vector of `PhyloTree` objects, to allow for
 multiple phylogenies, as will be the case when recombination is included.
 If not including recombination, `trees` will simply be of length 1.

 */


class PhyloOneChrom {

public:
    std::vector<PhyloTree> trees;
    std::vector<VarChrom*> var_chrom_ptrs;  // pointers to final VarChrom objects
    std::vector<VarChrom> tmp_chroms;       // temporary VarChrom's to evolve across tree
    std::vector<std::deque<uint8>> rates;   // rate indices (Gammas + invariants) for tree
    TreeMutator mutator;                    // to do the mutation additions across tree
    uint64 n_tips;                          // number of tips (i.e., variants)


    PhyloOneChrom() {}
    /*
     Construct just tree and mutator info, when no VarSet info yet available.
     Used in `fill_tree_mutator` method below.
     */
    PhyloOneChrom(
        const std::vector<uint64>& n_bases_,
        const std::vector<std::vector<double>>& branch_lens_,
        const std::vector<arma::Mat<uint64>>& edges_,
        const std::vector<std::vector<std::string>>& tip_labels_,
        const TreeMutator& mutator_base
    )
        : trees(edges_.size()),
          var_chrom_ptrs(),
          tmp_chroms(),
          mutator(mutator_base),
          chrom_rates(edges_[0].max()),
          n_tips(tip_labels_[0].size()),
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
        uint64 end_ = 0; // note: non-inclusive end point
        for (uint64 i = 0; i < n_trees; i++) {
            end_ += n_bases_[i];
            if (i > 0) start_ += n_bases_[i-1];
            trees[i] = PhyloTree(branch_lens_[i], edges_[i], tip_labels_[i],
                                 start_, end_);
        }

    }

    /*
     Set variant info:
     */
    void set_var_info(VarSet& var_set, const uint64& chrom_ind) {

        uint64 tree_size = trees[0].tree_size;
        uint64 n_vars = var_set.size();

        // Filling in pointers:
        var_chrom_ptrs.resize(n_vars);
        for (uint64 i = 0; i < n_vars; i++) {
            var_chrom_ptrs[i] = &var_set[i][chrom_ind];
        }

        // Fill in blank VarChrom objects:
        tmp_chroms = std::vector<VarChrom>(tree_size,
                                            VarChrom((*var_set.reference)[chrom_ind]));

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
     Fill tree and mutator info from an input list and base mutator object
    */
    void fill_tree_mutator(const List& genome_phylo_info, const uint64& i,
                           const TreeMutator& mutator_base) {

        std::string err_msg;

        const List& chrom_phylo_info(genome_phylo_info[i]);

        uint64 n_trees = chrom_phylo_info.size();
        if (n_trees == 0) {
            err_msg = "\nNo trees supplied on chromosome " + std::to_string(i+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        std::vector<uint64> n_bases_(n_trees);
        std::vector<std::vector<double>> branch_lens_(n_trees);
        std::vector<arma::Mat<uint64>> edges_(n_trees);
        std::vector<std::vector<std::string>> tip_labels_(n_trees);

        for (uint64 j = 0; j < n_trees; j++) {

            const List& phylo_info(chrom_phylo_info[j]);
            std::vector<double> branch_lens = as<std::vector<double>>(
                phylo_info["branch_lens"]);
            arma::Mat<uint64> edges = as<arma::Mat<uint64>>(phylo_info["edges"]);
            std::vector<std::string> tip_labels = as<std::vector<std::string>>(
                phylo_info["labels"]);
            if (branch_lens.size() != edges.n_rows) {
                err_msg = "\nBranch lengths and edges don't have the same ";
                err_msg += "size on chromosome " + std::to_string(i+1) + " and tree ";
                err_msg += std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }
            if (branch_lens.size() == 0) {
                err_msg = "\nEmpty tree on chromosome " + std::to_string(i+1);
                err_msg += " and tree " + std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }
            uint64 start = as<uint64>(phylo_info["start"]);
            uint64 end = as<uint64>(phylo_info["end"]);
            end++; // note: non-inclusive end point
            if (end <= start) {
                err_msg = "\nEnd position <= start position on chromosome ";
                err_msg += std::to_string(i+1) + " and tree " + std::to_string(j+1);
                throw(Rcpp::exception(err_msg.c_str(), false));
            }

            n_bases_[j] = end - start;
            branch_lens_[j] = branch_lens;
            edges_[j] = edges;
            tip_labels_[j] = tip_labels;
        }


        *this = PhyloOneChrom(n_bases_, branch_lens_, edges_, tip_labels_, mutator_base);

        return;

    }




private:


    bool recombination;


    /*
     Evolve one tree.
     */
    int one_tree(PhyloTree& tree, pcg64& eng, Progress& prog_bar);


    /*
     Reset for a new tree:
     */
    void reset(const PhyloTree& tree, pcg64& eng) {

        const uint64& tree_size(tree.tree_size);
        const uint64& start(tree.start);
        const uint64& end(tree.end);

        if (tree_size == 0) {
            throw(Rcpp::exception("\ntree size of zero is non-sensical.", false));
        }
        // Resize blank VarChrom objects if necessary:
        if (tree_size != tmp_chroms.size()) {
            VarChrom var_chrom_(*(var_chrom_ptrs[0]->ref_chrom));
            tmp_chroms.resize(tree_size, var_chrom_);
        }
        // Empty mutations from tree of VarChrom objects:
        for (VarChrom& var_chrom : tmp_chroms) var_chrom.clear();

        return;
    }


    /*
     Update for a new edge:
     */
    void update(std::exponential_distribution<double>& distr,
                const uint64& b1,
                const uint64& b2) {

        /*
         Replace existing mutation information in VarChrom at `b1` with info in the
         one at `b2`
         */
        samplers[b2].var_chrom->replace(*samplers[b1].var_chrom);
        /*
         Do the same for the ChromGammas in the sampler:
         */
        samplers[b2].location.regions = samplers[b1].location.regions;


        // Set exponential distribution to use this chromosome's rate:
        distr.param(std::exponential_distribution<double>::param_type(chrom_rates[b2]));

        return;
    }

    /*
     Clear info from VarChrom object at `b1` if it's no longer needed, to free up
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
        if (clear_b1) samplers[b1].var_chrom->clear();
        return;
    }

};







/*
 Phylogenetic tree info for all chromosomes in a genome.

 */
class PhyloInfo {
public:

    std::vector<PhyloOneChrom> phylo_one_chroms;

    PhyloInfo(const List& genome_phylo_info, const TreeMutator& mutator_base) {

        uint64 n_chroms = genome_phylo_info.size();

        if (n_chroms == 0) {
            throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                                  false));
        }

        phylo_one_chroms = std::vector<PhyloOneChrom>(n_chroms);

        // Fill tree and mutator info (i.e., everything but variant info):
        for (uint64 i = 0; i < n_chroms; i++) {
            phylo_one_chroms[i].fill_tree_mutator(genome_phylo_info, i, mutator_base);
        }
    }

    XPtr<VarSet> evolve_chroms(
            SEXP& ref_genome_ptr,
            const uint64& n_threads,
            const bool& show_progress);



};








#endif
