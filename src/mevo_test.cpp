
/*
 ****************************************************

 Functions for testing molecular evolution code.
 These functions are R versions of what would normally be run entirely from C++.
 This is to test the output in R using the testthat package.

 ****************************************************
 */


#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <progress.hpp>  // for the progress bar



#include "gemino_types.h"
#include "mevo.h"
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class
#include "mevo_phylo.h"  // match_ and template functions
#include "mevo_rate_matrices.h"  // rate matrix functions

using namespace Rcpp;




// Turn a Mutation into a List
List conv_mut(const Mutation& mut) {
    List out = List::create(_["size_modifier"] = mut.size_modifier,
                            _["old_pos"] = mut.old_pos,
                            _["new_pos"] = mut.new_pos,
                            _["nucleos"] = mut.nucleos);
    return out;
}



//' Turns a VarGenome's mutations into a list of data frames.
//'
//' Internal function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
DataFrame see_mutations(SEXP var_set_, const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_);
    const VarGenome& var_genome((*var_set)[var_ind]);

    uint32 n_muts = 0;
    for (const VarSequence& vs : var_genome.var_genome) n_muts += vs.mutations.size();

    std::vector<sint32> size_mod;
    size_mod.reserve(n_muts);
    std::vector<uint32> old_pos;
    old_pos.reserve(n_muts);
    std::vector<uint32> new_pos;
    new_pos.reserve(n_muts);
    std::vector<std::string> nucleos;
    nucleos.reserve(n_muts);
    std::vector<uint32> vars(n_muts, var_ind);
    std::vector<uint32> seqs;
    seqs.reserve(n_muts);

    for (uint32 i = 0; i < var_genome.size(); i++) {
        const VarSequence& var_seq(var_genome.var_genome[i]);
        uint32 n_muts_i = var_seq.mutations.size();
        for (uint32 j = 0; j < n_muts_i; ++j) {
            size_mod.push_back(var_seq.mutations[j].size_modifier);
            old_pos.push_back(var_seq.mutations[j].old_pos);
            new_pos.push_back(var_seq.mutations[j].new_pos);
            nucleos.push_back(var_seq.mutations[j].nucleos);
            seqs.push_back(i);
        }
    }

    DataFrame out = DataFrame::create(
        _["var"] = vars,
        _["seq"] = seqs,
        _["size_mod"] = size_mod,
        _["old_pos"] = old_pos,
        _["new_pos"] = new_pos,
        _["nucleos"] = nucleos);

    return out;
}


//' Turns a VarGenome's mutations into a list of data frames.
//'
//' Internal function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
List examine_mutations(SEXP var_set_, const uint32& var_ind, const uint32& seq_ind) {

    XPtr<VarSet> var_set_xptr(var_set_);
    const VarGenome& var_genome((*var_set_xptr)[var_ind]);
    const VarSequence& var_seq(var_genome[seq_ind]);

    std::string bases = "TCAG";
    std::vector<uint32> base_inds(85);
    uint32 j = 0;
    for (const char& c : bases) {
        base_inds[static_cast<uint32>(c)] = j;
        j++;
    }

    uint32 n_muts = var_seq.mutations.size();
    arma::mat sub_mat(4, 4, arma::fill::zeros);
    uint32 max_ins = 0;
    uint32 max_del = 0;
    for (uint32 i = 0; i < n_muts; i++) {
        sint32 mi = var_seq.mutations[i].size_modifier;
        if (mi == 0) continue;
        if (mi > 0) {
            if (mi > max_ins) max_ins = mi;
        } else {
            uint32 mid = static_cast<uint32>(std::abs(mi));
            if (mid > max_del) max_del = mid;
        }
    }
    arma::mat ins_mat(4, max_ins, arma::fill::zeros);
    arma::mat del_mat(4, max_del, arma::fill::zeros);
    std::vector<uint32> pos_vec(n_muts);

    for (uint32 mut_i = 0; mut_i < n_muts; mut_i++) {

        const Mutation& m(var_seq.mutations[mut_i]);

        char c = var_seq.ref_seq[m.old_pos];
        uint32 i = base_inds[static_cast<uint32>(c)];
        sint32 smod = m.size_modifier;
        if (smod == 0) {
            uint32 j = base_inds[static_cast<uint32>(m.nucleos[0])];
            sub_mat(i, j)++;
        } else if (smod > 0) {
            uint32 j = static_cast<uint32>(smod - 1);
            ins_mat(i, j)++;
        } else {
            uint32 j = static_cast<uint32>(std::abs(smod + 1));
            del_mat(i, j)++;
        }

        pos_vec[mut_i] = var_seq.mutations[mut_i].old_pos;
    }

    List out = List::create(
        _["sub"] = wrap(sub_mat),
        _["ins"] = wrap(ins_mat),
        _["del"] = wrap(del_mat),
        _["pos"] = pos_vec);

    return out;
}


//' Faster version of table function to count the number of mutations in Gamma regions.
//'
//' @param gamma_ends Vector of endpoints for gamma regions
//' @param positions Vector of positions that you want to bin into gamma regions.
//'
//'
//'
//[[Rcpp::export]]
std::vector<uint32> table_gammas(const std::vector<uint32>& gamma_ends,
                                 const std::vector<uint32>& positions) {
    std::vector<uint32> out(gamma_ends.size(), 0U);
    for (uint32 i = 0; i < positions.size(); i++) {
        uint32 j = std::lower_bound(gamma_ends.begin(), gamma_ends.end(),
                                  positions[i]) - gamma_ends.begin();
        out[j]++;
    }
    return out;
}



//' Add mutations manually from R.
//'
//' Note that all indices are in 0-based C++ indexing. This means that the first
//' item is indexed by `0`, and so forth.
//'
//' @param var_set_ External pointer to a C++ `VarSet` object
//' @param var_ind Integer index to the desired variant. Uses 0-based indexing!
//' @param seq_ind Integer index to the desired sequence. Uses 0-based indexing!
//' @param new_pos_ Integer index to the desired subsitution location.
//'     Uses 0-based indexing!
//'
//' @name add_mutations
NULL_ENTRY;

//' @describeIn add_mutations Add a substitution.
//'
//' @inheritParams var_set_ add_mutations
//' @inheritParams var_ind add_mutations
//' @inheritParams seq_ind add_mutations
//' @param nucleo_ Character to substitute for existing one.
//' @inheritParams new_pos_ add_mutations
//'
//'
//[[Rcpp::export]]
void add_substitution(SEXP var_set_, const uint32& var_ind,
                      const uint32& seq_ind,
                      const char& nucleo_,
                      const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_);
    VarGenome& var_genome((*var_set)[var_ind]);
    VarSequence& var_seq(var_genome[seq_ind]);
    var_seq.add_substitution(nucleo_, new_pos_);
    return;
}
//' @describeIn add_mutations Add an insertion.
//'
//' @inheritParams var_set_ add_mutations
//' @inheritParams var_ind add_mutations
//' @inheritParams seq_ind add_mutations
//' @param nucleos_ Nucleotides to insert at the desired location.
//' @inheritParams new_pos_ add_mutations
//'
//'
//[[Rcpp::export]]
void add_insertion(SEXP var_set_, const uint32& var_ind,
                   const uint32& seq_ind,
                   const std::string& nucleos_,
                   const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_);
    VarGenome& var_genome((*var_set)[var_ind]);
    VarSequence& var_seq(var_genome[seq_ind]);
    var_seq.add_insertion(nucleos_, new_pos_);
    return;
}
//' @describeIn add_mutations Add a deletion.
//'
//' @inheritParams var_set_ add_mutations
//' @inheritParams var_ind add_mutations
//' @inheritParams seq_ind add_mutations
//' @param size_ Size of deletion.
//' @inheritParams new_pos_ add_mutations
//'
//'
//[[Rcpp::export]]
void add_deletion(SEXP var_set_, const uint32& var_ind,
                  const uint32& seq_ind,
                  const uint32& size_,
                  const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_);
    VarGenome& var_genome((*var_set)[var_ind]);
    VarSequence& var_seq(var_genome[seq_ind]);
    var_seq.add_deletion(size_, new_pos_);
    return;
}




//' Get a rate for given start and end points of a VarSequence.
//'
//' @noRd
//'
//[[Rcpp::export]]
double test_rate(const uint32& start, const uint32& end,
                 const uint32& var_ind, const uint32& seq_ind,
                 SEXP var_set_, SEXP sampler_,
                 const arma::mat& gamma_mat_) {

    XPtr<VarSet> var_set(var_set_);
    VarSequence& var_seq((*var_set)[var_ind][seq_ind]);

    XPtr<ChunkMutationSampler> sampler(sampler_);

    sampler->fill_ptrs(var_seq);
    sampler->fill_gamma(gamma_mat_);

    double out = sampler->total_rate(start, end, true);

    return out;

}






//' Test sampling based on an evolutionary model.
//'
//' Make SURE `sampler_base_` is a `ChunkMutationSampler`, not a `MutationSampler`!
//'
//' @param var_set_ Pointer to a VarSet object.
//' @param sampler_base_ Pointer to a ChunkMutationSampler object.
//' @param branch_lens Branch lengths from phylogeny.
//' @param edges Edge matrix from phylogeny.
//' @param tip_labels Character vector of the actual phylogeny's tip labels.
//' @param ordered_tip_labels Character vector of the tip labels in the order
//'     you want them.
//' @param gamma_mat Gamma matrix.
//' @param recombination Boolean for whether to include recombination. If this is
//'     \code{FALSE}, then \code{start} and \code{end} arguments are ignored.
//'     Defaults to \code{FALSE}.
//' @param start Starting point of region in which to insert mutations.
//'     Ignored if \code{recombination} is \code{FALSE}.
//' @param end Ending point of region in which to insert mutations.
//'     Ignored if \code{recombination} is \code{FALSE}.
//'
//'
//' @return A vector of integers indicating the number of mutations per edge.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::vector<uint32>> test_mevo(
        SEXP& var_set_,
        SEXP& sampler_base_,
        const std::vector<uint32>& seq_inds,
        const std::vector<double>& branch_lens,
        arma::Mat<uint32> edges,
        const std::vector<std::string>& tip_labels,
        const std::vector<std::string>& ordered_tip_labels,
        const std::vector<arma::mat>& gamma_mats,
        const bool& recombination = false,
        const uint32& start = 0,
        const sint64& end = 0) {

    XPtr<VarSet> var_set(var_set_);
    XPtr<ChunkMutationSampler> sampler_base(sampler_base_);

    uint32 n_seqs = seq_inds.size();

    if (n_seqs != gamma_mats.size()) {
        stop("seq_inds and gamma_mats must be the same length");
    }

    uint32 n_tips = var_set->size();
    if (ordered_tip_labels.size() != n_tips || tip_labels.size() != n_tips) {
        stop("ordered_tip_labels and tip_labels must have lengths == # variants.");
    }

    std::vector<uint32> spp_order = match_(ordered_tip_labels, tip_labels);

    uint32 n_edges = edges.n_rows;
    if (branch_lens.size() != n_edges) {
        stop("branch_lens must have the same length as the # rows in edges.");
    }
    if (edges.n_cols != 2) stop("edges must have exactly two columns.");
    // From R to C++ indices
    edges -= 1;

    std::vector<std::vector<uint32>> out(n_seqs, std::vector<uint32>(n_edges, 0));

    pcg32 eng = seeded_pcg();

    for (uint32 i = 0; i < n_seqs; i++) {
        std::vector<uint32>& n_muts(out[i]);
        const arma::mat& gamma_mat(gamma_mats[i]);
        const uint32& seq_ind(seq_inds[i]);
        if (gamma_mat(gamma_mat.n_rows-1,0) != (*var_set)[0][seq_ind].size()) {
            stop("gamma_mat doesn't reach the end of the sequence.");
        }

        int code = one_tree_<ChunkMutationSampler>(
            *var_set, *sampler_base, seq_ind, branch_lens, edges, spp_order,
            gamma_mat, eng, n_muts, recombination, start, end
        );

        // Make sure this happens outside of multithreaded code
        if (code == -1) {
            std::string warn_msg = "\nUser interrupted phylogenetic evolution. ";
            warn_msg += "Note that changes occur in place, so your variants have ";
            warn_msg += "already been partially added.";
            // throw(Rcpp::exception(err_msg.c_str(), false));
            Rcpp::warning(warn_msg.c_str());
            return out;
        }
    }

    return out;
}







