
/*
 ********************************************************

 Retrieving and printing info from Ref* and Var* classes, from R environment.
 Many of these functions are used for testing.

 ********************************************************
 */


#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng


#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "mevo.h"
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class
#include "mevo_phylo.h"  // match_ and template functions



using namespace Rcpp;



/*
 ========================================================================================
 ========================================================================================

 Printing

 ========================================================================================
 ========================================================================================
 */



//' Function to print info on a `RefGenome`.
//'
//' Access `RefGenome` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_ref_genome(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    ref_genome->print();
    return;
}



//' Function to print info on a VarSet.
//'
//' Access `VarSet` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_var_set(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    var_set->print();
    return;
}




/*
 ========================================================================================
 ========================================================================================

 Making new

 ========================================================================================
 ========================================================================================
 */



//' Make a RefGenome object from a set of sequences.
//'
//' Used for testing.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_ref_genome(const std::deque<std::string>& seqs) {
    XPtr<RefGenome> ref_genome(new RefGenome(seqs), true);
    return ref_genome;
}


//' Make a VarSet object from a RefGenome pointer and # variants.
//'
//' Used for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_var_set(SEXP ref_genome_ptr, const uint32& n_vars) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    XPtr<VarSet> var_set(new VarSet(*ref_genome, n_vars), true);
    return var_set;
}




/*
 ========================================================================================
 ========================================================================================

 Viewing numbers of sequences/variants

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
uint32 view_ref_genome_nseqs(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    uint32 out = ref_genome->size();
    return out;
}


// Number of sequences
//[[Rcpp::export]]
uint32 view_var_set_nseqs(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    uint32 out = var_set->reference->size();
    return out;
}

// Number of variants
//[[Rcpp::export]]
uint32 view_var_set_nvars(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    uint32 out = var_set->size();
    return out;
}


/*
 ========================================================================================
 ========================================================================================

 Viewing sequence sizes

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
std::vector<uint32> view_ref_genome_seq_sizes(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<uint32> out = ref_genome->seq_sizes();
    return out;
}


//' See all sequence sizes in a VarGenome object within a VarSet.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<uint32> view_var_genome_seq_sizes(SEXP var_set_ptr,
                                             const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
    const VarGenome& var_genome((*var_set)[var_ind]);

    std::vector<uint32> out(var_genome.size());
    for (uint32 i = 0; i < var_genome.size(); i++) {
        const VarSequence& var_seq(var_genome.var_genome[i]);
        out[i] = var_seq.seq_size;
    }
    return out;
}







/*
 ========================================================================================
 ========================================================================================

 Viewing one sequence from a genome

 ========================================================================================
 ========================================================================================
 */



//[[Rcpp::export]]
std::string view_ref_genome_seq(SEXP ref_genome_ptr, const uint32& seq_ind) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::string out = (*ref_genome)[seq_ind].nucleos;
    return out;
}


//' Function to piece together the strings for one sequence in a VarGenome.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::string view_var_genome_seq(SEXP var_set_ptr,
                               const uint32& var_ind,
                               const uint32& seq_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
    const VarSequence& var_seq((*var_set)[var_ind][seq_ind]);
    std::string out = var_seq.get_seq_full();
    return out;
}




/*
 ========================================================================================
 ========================================================================================

 Viewing an entire genome

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
std::vector<std::string> view_ref_genome(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<std::string> out(ref_genome->size(), "");
    for (uint32 i = 0; i < ref_genome->size(); i++) {
        const RefSequence& ref_seq((*ref_genome)[i]);
        out[i] = ref_seq.nucleos;
    }
    return out;
}

//' Function to piece together the strings for all sequences in a VarGenome.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> view_var_genome(SEXP var_set_ptr,
                                        const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
    const VarGenome& var_genome((*var_set)[var_ind]);

    std::vector<std::string> out(var_genome.size(), "");
    for (uint32 i = 0; i < var_genome.size(); i++) {
        const VarSequence& var_seq(var_genome[i]);
        out[i] = var_seq.get_seq_full();
    }
    return out;
}





/*
 ========================================================================================
 ========================================================================================

 Viewing names

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
std::vector<std::string> view_ref_genome_seq_names(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<std::string> out;
    out.reserve(ref_genome->size());
    for (const RefSequence& seq : (*ref_genome).sequences) out.push_back(seq.name);
    return out;
}

//' See all variant-genome names in a VarSet object.
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<std::string> view_var_set_var_names(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    std::vector<std::string> out;
    out.reserve(var_set->size());
    for (const VarGenome& vg : var_set->variants) out.push_back(vg.name);
    return out;
}




/*
 ========================================================================================
 ========================================================================================

 Setting names

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
void set_ref_genome_seq_names(
        SEXP ref_genome_ptr,
        const std::vector<uint32>& seq_inds,
        const std::vector<std::string>& names) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    if (names.size() != seq_inds.size()) stop("names and seq_inds aren't the same size");
    if (*std::max_element(seq_inds.begin(), seq_inds.end()) >= ref_genome->size()) {
        stop("at least one value in seq_inds is too large");
    }
    for (uint32 i = 0; i < seq_inds.size(); i++) {
        (*ref_genome)[seq_inds[i]].name = names[i];
    }
    return;
}


//[[Rcpp::export]]
void set_var_set_var_names(
        SEXP var_set_ptr,
        const std::vector<uint32>& var_inds,
        const std::vector<std::string>& names) {
    XPtr<VarSet> var_set(var_set_ptr);
    if (names.size() != var_inds.size()) stop("names and var_inds aren't the same size");
    if (*std::max_element(var_inds.begin(), var_inds.end()) >= var_set->size()) {
        stop("at least one value in var_inds is too large");
    }
    for (uint32 i = 0; i < var_inds.size(); i++) {
        (*var_set)[var_inds[i]].name = names[i];
    }
    return;
}



/*
 ========================================================================================
 ========================================================================================

 Remove elements

 ========================================================================================
 ========================================================================================
 */


//[[Rcpp::export]]
void remove_ref_genome_seqs(
        SEXP ref_genome_ptr,
        std::vector<uint32> seq_inds) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefSequence>& sequences(ref_genome->sequences);

    // Checking for duplicates:
    std::sort(seq_inds.begin(), seq_inds.end());
    if (adjacent_find(seq_inds.begin(), seq_inds.end()) != seq_inds.end()) {
        stop("duplicates detected in seq_inds");
    }

    // Number of deleted nucleotides:
    uint64 n_del = 0;

    for (uint32 i = 1; i <= seq_inds.size(); i++) {
        // Going backward so I don't have to update later ones each time:
        uint32 j = seq_inds[(seq_inds.size() - i)];
        n_del += sequences[j].size();
        sequences.erase(sequences.begin() + j);
    }
    clear_memory<std::deque<RefSequence>>(sequences);

    // Update total size:
    ref_genome->total_size -= n_del;

    return;
}



//[[Rcpp::export]]
void remove_var_set_vars(
        SEXP var_set_ptr,
        std::vector<uint32> var_inds) {

    XPtr<VarSet> var_set(var_set_ptr);
    std::deque<VarGenome>& variants(var_set->variants);

    // Checking for duplicates:
    std::sort(var_inds.begin(), var_inds.end());
    if (adjacent_find(var_inds.begin(), var_inds.end()) != var_inds.end()) {
        stop("duplicates detected in var_inds");
    }

    for (uint32 i = 1; i <= var_inds.size(); i++) {
        // Going backward so I don't have to update later ones each time:
        uint32 j = var_inds[(var_inds.size() - i)];
        variants.erase(variants.begin() + j);
    }
    clear_memory<std::deque<VarGenome>>(variants);
    return;
}







/*
 ========================================================================================
 ========================================================================================

 For molecular evolution testing

 ========================================================================================
 ========================================================================================
 */


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
DataFrame view_mutations(SEXP var_set_ptr, const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
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
List examine_mutations(SEXP var_set_ptr, const uint32& var_ind, const uint32& seq_ind) {

    XPtr<VarSet> var_set_xptr(var_set_ptr);
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

        char c = (*(var_seq.ref_seq))[m.old_pos];
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
//' @param var_set_ptr External pointer to a C++ `VarSet` object
//' @param var_ind Integer index to the desired variant. Uses 0-based indexing!
//' @param seq_ind Integer index to the desired sequence. Uses 0-based indexing!
//' @param new_pos_ Integer index to the desired subsitution location.
//'     Uses 0-based indexing!
//'
//' @name add_mutations
NULL_ENTRY;

//' @describeIn add_mutations Add a substitution.
//'
//' @inheritParams add_mutations
//' @param nucleo_ Character to substitute for existing one.
//'
//'
//[[Rcpp::export]]
void add_substitution(SEXP var_set_ptr, const uint32& var_ind,
                      const uint32& seq_ind,
                      const char& nucleo_,
                      const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_ptr);
    VarGenome& var_genome((*var_set)[var_ind]);
    VarSequence& var_seq(var_genome[seq_ind]);
    var_seq.add_substitution(nucleo_, new_pos_);
    return;
}
//' @describeIn add_mutations Add an insertion.
//'
//' @inheritParams add_mutations
//' @param nucleos_ Nucleotides to insert at the desired location.
//'
//'
//[[Rcpp::export]]
void add_insertion(SEXP var_set_ptr, const uint32& var_ind,
                   const uint32& seq_ind,
                   const std::string& nucleos_,
                   const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_ptr);
    VarGenome& var_genome((*var_set)[var_ind]);
    VarSequence& var_seq(var_genome[seq_ind]);
    var_seq.add_insertion(nucleos_, new_pos_);
    return;
}
//' @describeIn add_mutations Add a deletion.
//'
//' @inheritParams add_mutations
//' @param size_ Size of deletion.
//'
//'
//[[Rcpp::export]]
void add_deletion(SEXP var_set_ptr, const uint32& var_ind,
                  const uint32& seq_ind,
                  const uint32& size_,
                  const uint32& new_pos_) {
    XPtr<VarSet> var_set(var_set_ptr);
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
                 SEXP var_set_ptr, SEXP sampler_base_ptr,
                 const arma::mat& gamma_mat_) {

    XPtr<VarSet> var_set(var_set_ptr);

    VarSequence& var_seq((*var_set)[var_ind][seq_ind]);

    XPtr<ChunkMutationSampler> sampler_base(sampler_base_ptr);

    ChunkMutationSampler sampler(*sampler_base);
    sampler.fill_ptrs(var_seq);
    sampler.fill_gamma(gamma_mat_);

    double out = sampler.total_rate(start, end, true);

    return out;

}


