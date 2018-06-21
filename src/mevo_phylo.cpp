
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
inline int one_tree_no_recomb_(VarSet& vars,
                               const T& sampler_base,
                               const uint32& seq_ind,
                               const std::vector<double>& branch_lens,
                               const arma::Mat<uint32>& edges,
                               const std::vector<uint32>& spp_order,
                               const arma::mat& gamma_mat,
                               pcg32& eng,
                               Progress& progress,
                               const std::vector<uint32>& progress_branch_lens,
                               std::vector<uint32>& n_muts) {


    // # tips = # variants
    uint32 n_tips = vars.size();
    // tree size is # tips plus # nodes
    uint32 tree_size = edges.max() + 1;
    // Number of edges = the number of connections between nodes/tips
    uint32 n_edges = edges.n_rows;

    /*
     Create tree of empty VarSequence objects
     */
    std::vector<VarSequence> var_seqs(tree_size, VarSequence(vars.reference[seq_ind]));

    /*
     Create corresponding tree of MutationSampler objects
     First fill with the base sampler, then fill in pointers to the corresponding
     VarSequence object, then fill in the gamma matrix.
     */
    std::vector<T> samplers(tree_size, sampler_base);
    for (uint32 i = 0; i < tree_size; i++) {
        samplers[i].fill_ptrs(var_seqs[i]);
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
    for (uint32 i = 0; i < n_edges; i++) {

        if (progress.is_aborted()) return -1;

        uint32 n_muts_ = 0;

        // Indices for nodes/tips that the branch length in `branch_lens` refers to
        uint32 b1 = edges(i,0);
        uint32 b2 = edges(i,1);
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
        double& rate(seq_rates[b2]);
        rate = seq_rates[b1];

        // Set exponential distribution to use this sequence's rate:
        distr.param(std::exponential_distribution<double>::param_type(rate));

        /*
         Now do exponential jumps and mutate until you exceed the branch length.
         */
        double amt_time = branch_lens[i];
        double time_jumped = distr(eng);
        double rate_change = 0;
        while (time_jumped <= amt_time && var_seqs[b2].size() > 0) {
            /*
             Add mutation here, outputting how much the overall sequence rate should
             change:
             */
            rate_change = samplers[b2].mutate_rate_change(eng);
            n_muts_++;
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
        if (i < (n_edges - 1)) {
            // Is it absent from any remaining items in the first column?
            clear_b1 = ! arma::any(edges(arma::span(i+1, edges.n_rows - 1), 0) == b1);
        } else clear_b1 = true;
        if (clear_b1) samplers[b1].vs->clear();

        progress.increment(progress_branch_lens[i]);

        n_muts[i] = n_muts_;
    }

    /*
     Update final `VarSequence` objects using the index vector `spp_order`:
     */
    for (uint32 i = 0; i < n_tips; i++) {
        uint32 j = spp_order[i];
        vars[i][seq_ind].replace(var_seqs[j]);
    }

    return 0;

}




/*
 This function produces `spp_order` in `test_phylo`.
 This object is a vector of which indices in the phylogeny tips should go first.
 This effectively ensures that the tip labels line up with the output from this function.
 It's equivalent to the following R code, where `ordered_tip_labels` is a character vector
 of the tip names in the order you always want them:
 `spp_order <- match(ordered_tip_labels, phy$tip.label)`
 */
std::vector<uint32> match_(std::vector<std::string> x, std::vector<std::string> table) {

    std::vector<uint32> out(x.size());

    for (uint32 i = 0; i < out.size(); i++) {
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
//' @return A vector of integers indicating the number of mutations per edge.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<uint32> test_phylo(SEXP& vs_sexp,
                               SEXP& sampler_base_sexp,
                               const uint32& seq_ind,
                               const std::vector<double>& branch_lens,
                               arma::Mat<uint32> edges,
                               const std::vector<std::string>& tip_labels,
                               const std::vector<std::string>& ordered_tip_labels,
                               const arma::mat& gamma_mat,
                               const bool& display_progress = false) {

    XPtr<VarSet> vs_xptr(vs_sexp);
    XPtr<ChunkMutationSampler> sampler_base_xptr(sampler_base_sexp);

    /*
     Setting up a sum of branch lengths to use for the progress bar.
     Since it gets cast to an integer, I process each branch length to be an integer >= 1.
     */
    std::vector<uint32> progress_branch_lens(branch_lens.size());
    double branch_min = *std::min_element(branch_lens.begin(), branch_lens.end());
    for (uint32 i = 0; i < progress_branch_lens.size(); i++) {
        progress_branch_lens[i] = static_cast<uint32>(branch_lens[i] / branch_min);
    }
    uint32 total_branches = std::accumulate(progress_branch_lens.begin(),
                                          progress_branch_lens.end(), 0);

    Progress progress(total_branches, display_progress);

    if (gamma_mat(gamma_mat.n_rows-1,0) != (*vs_xptr)[0][seq_ind].size()) {
        stop("gamma_mat doesn't reach the end of the sequence.");
    }

    uint32 n_tips = vs_xptr->size();
    if (ordered_tip_labels.size() != n_tips || tip_labels.size() != n_tips) {
        stop("ordered_tip_labels and tip_labels must have the same length as ",
             "# variants.");
    }

    std::vector<uint32> spp_order = match_(ordered_tip_labels, tip_labels);

    uint32 n_edges = edges.n_rows;
    if (branch_lens.size() != n_edges) {
        stop("branch_lens must have the same length as the # rows in edges.");
    }
    if (edges.n_cols != 2) stop("edges must have exactly two columns.");

    std::vector<uint32> n_muts(n_edges, 0);

    // From R to C++ indices
    edges -= 1;

    pcg32 eng = seeded_pcg();

    int code = one_tree_no_recomb_<ChunkMutationSampler>(
        *vs_xptr, *sampler_base_xptr, seq_ind, branch_lens, edges, spp_order,
        gamma_mat, eng, progress, progress_branch_lens, n_muts
    );

    // Make sure this happens outside of multithreaded code
    if (code == -1) {
        std::string warn_msg = "\nUser interrupted phylogenetic evolution. ";
        warn_msg += "Note that changes occur in place, so your variants have ";
        warn_msg += "already been partially added.";
        // throw(Rcpp::exception(err_msg.c_str(), false));
        Rcpp::warning(warn_msg.c_str());
    }

    return n_muts;
}



//' Get a rate for given start and end points of a VarSequence.
//'
//' @noRd
//'
//[[Rcpp::export]]
double test_rate(const uint32& start, const uint32& end,
                 const uint32& var_ind, const uint32& seq_ind,
                 SEXP var_set_sexp, SEXP sampler_sexp) {

    XPtr<VarSet> var_set(var_set_sexp);
    VarSequence& vs((*var_set)[var_ind][seq_ind]);

    XPtr<ChunkMutationSampler> sampler(sampler_sexp);

    arma::mat gamma_mat(1, 2);
    gamma_mat(0,0) = vs.size();
    gamma_mat(0,1) = 1;

    sampler->fill_ptrs(vs);
    sampler->fill_gamma(gamma_mat);

    double out = sampler->total_rate(start, end);

    return out;

}



