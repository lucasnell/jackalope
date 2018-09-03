
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

using namespace Rcpp;



/*
 ------------------
 Retrieve a nucleotide (char type) from the variant sequence
 based on the position in the new, variant sequence
 ------------------
 */
char VarSequence::get_nt(const uint32& new_pos) const {
    char out;
    /*
     Index to the Mutation object nearest to (without being past)
     an input position (see below for the `get_mut_` fxn):
     */
    uint32 mut_i = get_mut_(new_pos);
    /*
     If the new_pos is less than the position for the first mutation
     or if mutations is empty
     (in which cases `get_mut_` returns `mutations.size()`),
     we just extract the character from the beginning of the reference string:
     */
    if (mut_i == mutations.size()) {
        out = (*ref_seq)[new_pos];
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

    if (mutations.empty()) return ref_seq->nucleos;

    // Index to the first Mutation object
    uint32 mut_i = 0;

    std::string out(seq_size, 'x');
    uint32 pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations[mut_i].new_pos) {
        out[pos] = (*ref_seq)[pos];
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint32 next_mut_i = mut_i + 1;
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
std::string VarSequence::get_seq_start(uint32 out_length) const {

    if (out_length > seq_size) out_length = seq_size;

    if (mutations.empty()) return ref_seq->nucleos.substr(0, out_length);

    std::string out(out_length, 'x');

    uint32 mut_i = 0;

    uint32 pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations[mut_i].new_pos) {
        out[pos] = (*ref_seq)[pos];
        if (pos == out_length - 1) return out;
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint32 next_mut_i = mut_i + 1;
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
                                const uint32& start,
                                const uint32& chunk_size,
                                uint32& mut_i) const {

    uint32 end = start + chunk_size - 1;

    if (start >= seq_size) {
        mut_i = mutations.size();
        chunk_str.clear();
        return;
    }
    // Making sure end doesn't go beyond the sequence bounds
    if (end >= seq_size) end = seq_size - 1;

    uint32 out_length = end - start + 1;

    // No need to mess around with mutations if there aren't any
    if (mutations.empty()) {
        chunk_str = ref_seq->nucleos.substr(start, out_length);
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

    uint32 pos = start;
    uint32 next_mut_i = mut_i + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `start` is before the first mutation)
     */
    while (pos < mutations[mut_i].new_pos && pos <= end) {
        chunk_str[pos - start] = (*ref_seq)[pos];
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
char VarSequence::get_char_(const uint32& new_pos,
                            const uint32& mut_i) const {
    const Mutation& m(mutations[mut_i]);
    char out;
    uint32 ind = new_pos - m.new_pos;
    if (static_cast<sint32>(ind) > m.size_modifier) {
        ind += (m.old_pos - m.size_modifier);
        out = (*ref_seq)[ind];
    } else {
        out = m[ind];
    }
    return out;
}






/*
 ========================================================================================
 ========================================================================================

 These functions are to access fields from R

 ========================================================================================
 ========================================================================================
 */


/*
 ----------------------------------------
 Reference genome methods
 ----------------------------------------
 */


//[[Rcpp::export]]
std::vector<uint32> see_ref_genome_seq_sizes(SEXP ref_genome_) {
    XPtr<RefGenome> ref_genome(ref_genome_);
    std::vector<uint32> out = ref_genome->seq_sizes();
    return out;
}
// Check for no duplicates
//[[Rcpp::export]]
void remove_ref_genome_seq_sizes(
        SEXP ref_genome_,
        std::vector<uint32> seq_inds) {

    XPtr<RefGenome> ref_genome(ref_genome_);
    std::deque<RefSequence>& sequences(ref_genome->sequences);

    std::sort(seq_inds.begin(), seq_inds.end());

    for (uint32 i = 0; i < seq_inds.size(); i++) {
        uint32 j = seq_inds[i];
        sequences.erase(sequences.begin() + j);
        for (uint32 k = i+1; k < seq_inds.size(); k++) seq_inds[k]--;
    }
    clear_memory<std::deque<RefSequence>>(sequences);
    return;
}


//[[Rcpp::export]]
std::vector<std::string> see_ref_genome_seq_names(SEXP ref_genome_) {
    XPtr<RefGenome> ref_genome(ref_genome_);
    std::vector<std::string> out;
    out.reserve(ref_genome->size());
    for (const RefSequence& seq : (*ref_genome).sequences) out.push_back(seq.name);
    return out;
}
// Check for seq_inds and names lengths being equal
//[[Rcpp::export]]
void set_ref_genome_seq_names(
        SEXP ref_genome_,
        const std::vector<uint32>& seq_inds,
        const std::vector<std::string>& names) {
    XPtr<RefGenome> ref_genome(ref_genome_);
    for (uint32 i = 0; i < seq_inds.size(); i++) {
        (*ref_genome)[seq_inds[i]].name = names[i];
    }
    return;
}

//[[Rcpp::export]]
std::string see_ref_genome_seq(SEXP ref_genome_, const uint32& seq_ind) {
    XPtr<RefGenome> ref_genome(ref_genome_);
    std::string out = (*ref_genome)[seq_ind].nucleos;
    return out;
}



/*
 ----------------------------------------
 Variants methods
 ----------------------------------------
 */



