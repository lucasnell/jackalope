
/*
 ********************************************************

 Basic methods for sequence classes: retrieving info and initializing.

 ********************************************************
 */


#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque


#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "read_write.h"     // reading fasta to VarSet

using namespace Rcpp;



/*
 ------------------
 Retrieve a nucleotide (char type) from the variant sequence
 based on the position in the new, variant sequence
 ------------------
 */
char VarSequence::get_nt(const uint& new_pos) const {
    char out;
    /*
     Index to the Mutation object nearest to (without being past)
     an input position (see below for the `get_mut_` fxn):
     */
    uint mut_i = get_mut_(new_pos);
    /*
     If the new_pos is less than the position for the first mutation
     or if mutations is empty
     (in which cases `get_mut_` returns `mutations.size()`),
     we just extract the character from the beginning of the reference string:
     */
    if (mut_i == mutations.size()) {
        out = ref_seq[new_pos];
        /*
         If not, then extract the character from the Mutation object that
         `mut` points to (see below for `get_char_` fxn):
         */
    } else {
        out = get_char_(new_pos, mut_i);
    }

    return out;
}



/*
 ------------------
 Retrieve all nucleotides (i.e., the full sequence; std::string type) from
 the variant sequence
 ------------------
 */

std::string VarSequence::get_seq_full() const {

    if (mutations.empty()) return ref_seq.nucleos;

    // Index to the first Mutation object
    uint mut_i = 0;

    std::string out(seq_size, 'x');
    uint pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations[mut_i].new_pos) {
        out[pos] = ref_seq[pos];
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint next_mut_i = mut_i + 1;
    while (next_mut_i < mutations.size()) {
        while (pos < mutations[next_mut_i].new_pos) {
            out[pos] = get_char_(pos, mut_i);
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos < seq_size) {
        out[pos] = get_char_(pos, mut_i);
        ++pos;
    }

    return out;
}




/*
 ------------------
 Retrieve the first part of a sequence from the variant sequence.
 ------------------
 */
std::string VarSequence::get_seq_start(uint out_length) const {

    if (out_length > seq_size) out_length = seq_size;

    if (mutations.empty()) return ref_seq.nucleos.substr(0, out_length);

    std::string out(out_length, 'x');

    uint mut_i = 0;

    uint pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations[mut_i].new_pos) {
        out[pos] = ref_seq[pos];
        if (pos == out_length - 1) return out;
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint next_mut_i = mut_i + 1;
    while (next_mut_i < mutations.size()) {
        while (pos < mutations[next_mut_i].new_pos) {
            out[pos] = get_char_(pos, mut_i);
            if (pos == out_length - 1) return out;
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos < seq_size) {
        out[pos] = get_char_(pos, mut_i);
        if (pos == out_length - 1) return out;
        ++pos;
    }

    return out;
}


/*
 ------------------
 Set an input string object to any chunk of a sequence from the variant sequence.
 Before anything, this function moves `mut` to the location right before this chunk's
 starting position. I keep this index around so I don't have to iterate through
 the entire mutation deque multiple times.
 If end position is beyond the size of the sequence, it changes `chunk_str` to the
 sequence from the start to the sequence end.
 If start position is beyond the size of the sequence, it sets `mut` to `mutations.end()`
 and clears `chunk_str`.
 ------------------
 */
void VarSequence::set_seq_chunk(std::string& chunk_str,
                                const uint& start,
                                const uint& chunk_size,
                                uint& mut_i) const {

    uint end = start + chunk_size - 1;

    if (start >= seq_size) {
        mut_i = mutations.size();
        chunk_str.clear();
        return;
    }
    // Making sure end doesn't go beyond the sequence bounds
    if (end >= seq_size) end = seq_size - 1;

    uint out_length = end - start + 1;

    // No need to mess around with mutations if there aren't any
    if (mutations.empty()) {
        chunk_str = ref_seq.nucleos.substr(start, out_length);
        return;
    }
    // Move mutation to the proper spot
    while (mut_i < mutations.size()) {
        if (start < mutations[mut_i].new_pos) break;
        ++mut_i;
    }
    if (mut_i != 0) --mut_i;
    // Adjust input string size if necessary
    if (chunk_str.size() != out_length) chunk_str.resize(out_length, 'x');

    uint pos = start;
    uint next_mut_i = mut_i + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `start` is before the first mutation)
     */
    while (pos < mutations[mut_i].new_pos && pos <= end) {
        chunk_str[pos - start] = ref_seq[pos];
        ++pos;
    }
    if (pos > end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `end`) but before the next mutation
     */
    while (next_mut_i < mutations.size()) {
        while (pos < mutations[next_mut_i].new_pos && pos <= end) {
            chunk_str[pos - start] = get_char_(pos, mut_i);
            ++pos;
        }
        if (pos > end) return;
        ++mut_i;
        ++next_mut_i;
    }

    /*
     If we reach the last mutation, add nucleotides until `end` (remember that above,
     I've made sure that `end < seq_size`).
     */
    while (pos <= end) {
        chunk_str[pos - start] = get_char_(pos, mut_i);
        ++pos;
    }

    return;
}


/*
 ------------------
 Internal function for finding character of either mutation or reference
 given an index (in the "new", variant sequence) and an index for a
 single Mutation object.
 This only works if you've already narrowed it down to the Mutation object
 that is directly previous to the index position.
 ------------------
 */
char VarSequence::get_char_(const uint& new_pos,
                            const uint& mut_i) const {
    const Mutation& m(mutations[mut_i]);
    char out;
    uint ind = new_pos - m.new_pos;
    if (static_cast<sint>(ind) > m.size_modifier) {
        ind += (m.old_pos - m.size_modifier);
        out = ref_seq[ind];
    } else {
        out = m[ind];
    }
    return out;
}






VarSet::VarSet(const std::string& fasta_file, const uint& n_vars,
               const bool& cut_names, const bool& remove_soft_mask) {
    RefGenome reference;
    fill_ref_noind(reference, fasta_file, cut_names, remove_soft_mask);

    VarGenome vg(reference);
    variants = std::deque<VarGenome>(n_vars, vg);
    for (uint i = 0; i < n_vars; i++) variants[i].name = "var" + std::to_string(i);

}
VarSet::VarSet(const std::string& fasta_file, const std::string& fai_file,
               const uint& n_vars,
               const bool& remove_soft_mask) {
    RefGenome reference;
    fill_ref_ind(reference, fasta_file, fai_file, remove_soft_mask);

    VarGenome vg(reference);
    variants = std::deque<VarGenome>(n_vars, vg);
    for (uint i = 0; i < n_vars; i++) variants[i].name = "var" + std::to_string(i);
}

VarSet::VarSet(const std::deque<std::string>& seqs, const uint& n_vars)
    : variants(), reference(seqs) {
    VarGenome vg(reference);
    variants = std::deque<VarGenome>(n_vars, vg);
    for (uint i = 0; i < n_vars; i++) variants[i].name = "var" + std::to_string(i);
}








/*
 ========================================================================================
 ========================================================================================
 ========================================================================================
 ========================================================================================

 Used for testing in `tests/testthat/test_mutations.R`.

 ========================================================================================
 ========================================================================================
 ========================================================================================
 ========================================================================================
 */


//' Make a VarSet object from a set of sequences and # variants
//'
//' @noRd
//[[Rcpp::export]]
SEXP make_vars(const std::deque<std::string>& seqs, const uint& n_vars) {
    XPtr<VarSet> vset(new VarSet(seqs, n_vars), true);
    return vset;
}

//' Function to piece together the strings for all sequences in a VarGenome.
//'
//' @noRd
//[[Rcpp::export]]
std::vector<std::string> see_vg(SEXP vs_, const uint& v) {

    XPtr<VarSet> vs(vs_);
    VarGenome& vg((*vs)[v]);

    std::vector<std::string> out(vg.size(), "");
    for (uint i = 0; i < vg.size(); i++) {
        const VarSequence& vs(vg[i]);
        std::string s = vs.get_seq_full();
        out[i] = s;
    }
    return out;
}

//' See all sequence sizes in a VarSet object.
//'
//' @noRd
//[[Rcpp::export]]
std::vector<uint> see_sizes(SEXP vs_, const uint& v) {

    XPtr<VarSet> vs(vs_);
    VarGenome& vg((*vs)[v]);

    std::vector<uint> out(vg.size());
    for (uint i = 0; i < vg.size(); i++) {
        const VarSequence& vs(vg.var_genome[i]);
        out[i] = vs.seq_size;
    }
    return out;
}






/*
 ========================================================================================
 ========================================================================================
 ========================================================================================
 ========================================================================================

 TEMPORARY FUNCTIONS FOR TESTING

 ========================================================================================
 ========================================================================================
 ========================================================================================
 ========================================================================================
 */




//[[Rcpp::export]]
SEXP make_ref(std::deque<std::string> input) {
    XPtr<RefGenome> ref(new RefGenome(input), true);
    return ref;
}

//[[Rcpp::export]]
std::string see_ref_seq(SEXP ref_, const uint& s) {
    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);
    std::string out(ref[s].nucleos);
    return out;
}

//' Get a reference genome sequence's name.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::string see_ref_name(SEXP ref_, const uint& s) {
    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);
    std::string out(ref[s].name);
    return out;
}

//' Get a reference genome sequence's size.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
uint see_ref_seq_size(SEXP ref_, const uint& s) {
    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);
    uint out = ref[s].nucleos.size();
    return out;
}

//' Get number of sequences in a reference genome.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
uint see_ref_n_seq(SEXP ref_) {
    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);
    uint out = ref.size();
    return out;
}

