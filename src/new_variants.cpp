
/*
 * THIS FILE SEEKS TO REPLACE THE OLD CLASSES
 */


#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque


#include "gemino_types.h"  // integer types
#include "new_variants.h"
#include "read_write.h"     // reading fasta to VarSet
#include "util.h"   // cpp_rando_seq
using namespace Rcpp;


/*
 Calling `base::options("width")$width`
*/
int get_width() {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function opt_r = base["options"];
    // Call the function and receive its list output
    List width_list = opt_r("width");
    int console_width = width_list["width"];
    return console_width;
}





// Printing reference genome info
void RefGenome::print() const {

    int console_width = get_width();

    // 32 characters is the narrowest I'll allow
    // (I'd start getting negative lengths and other problems otherwise)
    console_width = std::max(console_width, 32);

    int num_seqs = size();
    std::vector<int> inds;
    if (num_seqs <= 10) {
        for (int i = 0; i < num_seqs; i++) inds.push_back(i);
    } else {
        for (int i = 0; i < 5; i++) inds.push_back(i);
        inds.push_back(-1);
        for (int i = (num_seqs - 5 + 1); i < num_seqs; i++) inds.push_back(i);
    }

    Rcout.imbue(std::locale(""));
    Rcout << "< Set of " << num_seqs << " sequences >" << std::endl;
    Rcout << "# Total size: " << total_size << " bp" << std::endl;
    // Rcout << "# Sequences:" << std::endl;

    int ind_i, name_width = 10, length_width = 9;
    // Console width minus name width AND length width AND spaces between
    int seq_print_len = console_width - name_width - length_width - 2;
    // Number of chars print before and after elipses for a long string
    int before_elips = std::ceil((seq_print_len - 3) / 2);
    int after_elips = seq_print_len - 3 - before_elips;

    Rprintf("%-*s %s%-*s %*s\n", name_width, "  name",
            std::string(before_elips - 4, ' ').c_str(),
            seq_print_len - (before_elips - 4), "sequence",
            length_width, "length");

    for (int i = 0; i < inds.size(); i++) {
        ind_i = inds[i];
        if (ind_i == -1) {
            Rprintf("%-10s %-*s %9s\n", "...", seq_print_len, "...", "...");
            continue;
        }
        const RefSequence& rs(sequences[ind_i]);
        const std::string& name_i(rs.name);
        const std::string& seq_i(rs.nucleos);
        // Print name
        Rprintf("%-10.10s ", name_i.c_str());
        // Print sequence
        if (seq_i.size() > seq_print_len){
            for (int j = 0; j < before_elips; j++) Rcout << seq_i[j];
            Rcout << "...";
            for (int j = (seq_i.size() - after_elips); j < seq_i.size(); j++) {
                Rcout << seq_i[j];
            }
        } else {
            Rprintf("%-*s", seq_print_len, seq_i.c_str());
        }
        // Print width
        if (rs.size() > 999999999) {
            Rprintf(" %9.2E", rs.size());
        } else {
            Rprintf(" %9i", rs.size());
        }
        Rcout << std::endl;
    }
}



/*
 ------------------
 Re-calculate new positions (and total scaffold size)
 ------------------
 */
// For all Mutation objects
void VarSequence::calc_positions(bool sort_first) {

    if (sort_first) std::sort(mutations.begin(), mutations.end());

    sint modifier = 0;

    for (auto iter = mutations.begin(); iter != mutations.end(); ++iter) {
        (*iter).new_pos = (*iter).old_pos + modifier;
        modifier += (*iter).size_modifier;
    }
    // Updating full scaffold size
    seq_size = ref_seq.size() + modifier;

    return;
}
/*
 For all Mutation objects after a given Mutation object
 (this is for after you insert a NEW Mutation, where `iter` below points to
 that Mutation)
 */
void VarSequence::calc_positions(std::deque<Mutation>::iterator iter) {

    sint modifier = (*iter).size_modifier;
    ++iter;

    // Updating individual Mutation objects
    for (; iter != mutations.end(); ++iter) {
        (*iter).new_pos += modifier;
    }
    // Updating full scaffold size
    seq_size += modifier;

    return;
}
// Variant on above to use only indices
void VarSequence::calc_positions(uint i) {

    sint modifier = mutations[i].size_modifier;
    ++i;

    // Updating individual Mutation objects
    for (; i < mutations.size(); i++) {
        mutations[i].new_pos += modifier;
    }
    // Updating full scaffold size
    seq_size += modifier;

    return;
}
/*
 For all Mutation objects after AND INCLUDING a given Mutation object
 (this is for after you MERGE multiple Mutations, where `iter` below points to
 that merged Mutation and `modifier` refers to the net change in sequence size
 after the merge)
 */
void VarSequence::calc_positions(std::deque<Mutation>::iterator iter,
                                  const sint& modifier) {

    // Updating individual Mutation objects
    for (; iter != mutations.end(); ++iter) {
        (*iter).new_pos += modifier;
    }
    // Updating full scaffold size
    seq_size += modifier;

    return;
}

// Variant on above to use only indices
void VarSequence::calc_positions(uint i, const sint& modifier) {
    // Updating individual Mutation objects
    for (; i < mutations.size(); ++i) {
        mutations[i].new_pos += modifier;
    }
    // Updating full scaffold size
    seq_size += modifier;

    return;
}


/*
 ------------------
 Retrieve a nucleotide (char type) from the variant scaffold
 based on the position in the new, variant scaffold
 ------------------
 */
char VarSequence::get_nt(const uint& new_pos) const {
    char out;
    /*
     Iterator to the Mutation object nearest to (without being past)
     an input position (see below for the `get_mut_` fxn):
     */
    std::deque<Mutation>::iterator mut = get_mut_(new_pos);
    /*
     If the new_pos is less than the position for the first mutation
     (in which case `get_mut_` returns `mutations.end()`), we
     just extract the character from the beginning of the reference string:
     */
    if (mut == mutations.end()) {
        out = ref_seq[new_pos];
        /*
         If not, then extract the character from the Mutation object that
         `mut` points to (see below for `get_char_` fxn):
         */
    } else {
        out = get_char_(new_pos, mut);
    }

    return out;
}





/*
 ------------------
 Retrieve all nucleotides (i.e., the full sequence; std::string type) from
 the variant scaffold
 ------------------
 */

std::string VarSequence::get_seq_full() const {

    if (mutations.empty()) return ref_seq.nucleos;

    // Iterator to the first Mutation object
    std::deque<Mutation>::iterator mut = mutations.begin();

    std::string out(seq_size, 'x');
    uint pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < (*mut).new_pos) {
        out[pos] = ref_seq[pos];
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    std::deque<Mutation>::const_iterator next_mut = mut + 1;
    while (next_mut != mutations.end()) {
        while (pos < (*next_mut).new_pos) {
            out[pos] = get_char_(pos, mut);
            ++pos;
        }
        ++mut;
        ++next_mut;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos < seq_size) {
        out[pos] = get_char_(pos, mut);
        ++pos;
    }

    return out;
}




/*
 ------------------
 Retrieve the first part of a sequence from the variant scaffold.
 ------------------
 */
std::string VarSequence::get_seq_start(uint out_length) const {

    if (out_length > seq_size) out_length = seq_size;

    if (mutations.empty()) return ref_seq.nucleos.substr(0, out_length);

    std::string out(out_length, 'x');

    std::deque<Mutation>::iterator mut = mutations.begin();

    uint pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < (*mut).new_pos) {
        out[pos] = ref_seq[pos];
        if (pos == out_length - 1) return out;
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    std::deque<Mutation>::const_iterator next_mut = mut + 1;
    while (next_mut != mutations.end()) {
        while (pos < (*next_mut).new_pos) {
            out[pos] = get_char_(pos, mut);
            if (pos == out_length - 1) return out;
            ++pos;
        }
        ++mut;
        ++next_mut;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos < seq_size) {
        out[pos] = get_char_(pos, mut);
        if (pos == out_length - 1) return out;
        ++pos;
    }

    return out;
}


/*
 ------------------
 Set an input string object to any chunk of a sequence from the variant scaffold.
 Before anything, this function moves `mut` to the location right before this chunk's
 starting position. I keep this iterator around so I don't have to iterate through
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
                                std::deque<Mutation>::iterator& mut) const {

    uint end = start + chunk_size - 1;

    if (start >= seq_size) {
        mut = mutations.end();
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
    while (mut != mutations.end()) {
        if (start < (*mut).new_pos) break;
        ++mut;
    }
    if (mut != mutations.begin()) --mut;
    // Adjust input string size if necessary
    if (chunk_str.size() != out_length) chunk_str.resize(out_length, 'x');

    uint pos = start;
    std::deque<Mutation>::const_iterator next_mut = mut + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `start` is before the first mutation)
     */
    while (pos < (*mut).new_pos && pos <= end) {
        chunk_str[pos - start] = ref_seq[pos];
        ++pos;
    }
    if (pos > end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `end`) but before the next mutation
     */
    while (next_mut != mutations.end()) {
        while (pos < (*next_mut).new_pos && pos <= end) {
            chunk_str[pos - start] = get_char_(pos, mut);
            ++pos;
        }
        if (pos > end) return;
        ++mut;
        ++next_mut;
    }

    /*
     If we reach the last mutation, add nucleotides until `end` (remember that above,
     I've made sure that `end < seq_size`).
     */
    while (pos <= end) {
        chunk_str[pos - start] = get_char_(pos, mut);
        ++pos;
    }

    return;
}




/*
 ------------------
 Add a deletion somewhere in the deque
 ------------------
 */
void VarSequence::add_deletion(const uint& size_, const uint& new_pos_) {

    std::deque<Mutation>::iterator iter;

    // Renaming this for a more descriptive name and to allow it to change
    uint deletion_start = new_pos_;

    /*
     Last position this deletion refers to
     (`std::min` is to coerce this deletion to one that's possible):
     */
    uint deletion_end = std::min(deletion_start + size_ - 1, seq_size - 1);

    /*
     Position on the old reference sequence for this deletion.
     This will change if there are insertions or deletions before `new_pos_`.
     */
    uint old_pos_ = new_pos_;

    /*
     Size modifier of this deletion. This can change when deletion merges with an
     insertion or another deletion.
     */
    sint size_mod = deletion_start - deletion_end - 1;


    /*
     If the mutations deque isn't empty, we may need to edit some of the mutations
     after this deletion if they're affected by it.
     See `deletion_blowup_` below for more info.
     */
    if (!mutations.empty()) {

        /*
         Scaffold-size modifier to be used to edit subsequent mutations.
         This number does not change.
         */
        const sint subseq_modifier(size_mod);

        iter = get_mut_(deletion_start);

        /*
         This "blows up" subsequent mutations if they're destroyed/altered
         by this deletion.
         See `deletion_blowup_` below for more info.
         */
        deletion_blowup_(iter, deletion_start, deletion_end, size_mod);

        /*
         If `size_mod` is zero, this means that an insertion/insertions absorbed all
         of the deletion, so after adjusting sizes, our business is done here.
         */
        if (size_mod == 0) {
            calc_positions(iter, subseq_modifier);
            return;
        }

        /*
         If the deletion hasn't been absorbed, we need to calculate its
         position on the old (i.e., reference) sequence:
         */
        if (iter != mutations.begin()) {
            --iter;
            old_pos_ = deletion_start - (*iter).new_pos + (*iter).old_pos -
                (*iter).size_modifier;
            ++iter;
        } else old_pos_ = deletion_start; // (`deletion_start` may have changed)

        // Adjust (1) positions of all mutations after and including `iter`, and
        //        (2) the scaffold size
        calc_positions(iter, subseq_modifier);
    /*
     If `mutations` is empty, just point to the beginning and adjust scaffold size
     */
    } else {
        iter = mutations.begin();
        seq_size += size_mod;
    }
    /*
     Now create the Mutation and insert it.
     */
    Mutation new_mut(old_pos_, deletion_start, size_mod);
    mutations.insert(iter, new_mut);

    return;
}




/*
 ------------------
 Add an insertion somewhere in the deque
 ------------------
 */
void VarSequence::add_insertion(const std::string& nucleos_, const uint& new_pos_) {
    // std::deque<Mutation>::iterator iter = get_mut_(new_pos_);
    uint i = get_mut_ind_(new_pos_);
    // `mutations.end()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    // if (iter == mutations.end()) {
    if (i == mutations.size()) {
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nucleos_);
        mutations.push_front(new_mut);
        // Adjust new positions and total scaffold size:
        calc_positions(static_cast<uint>(0));
        return;
    }
    // uint ind = new_pos_ - (*iter).new_pos;
    uint ind = new_pos_ - mutations[i].new_pos;
    /*
     If `new_pos_` is within the Mutation sequence, we adjust that Mutation:
     */
    // if (ind <= (*iter).size_modifier) {
    if (ind <= mutations[i].size_modifier) {
        sint size_ = nucleos_.size() - 1;
        // string to store combined nucleotides
        std::string nt = "";
        // for (uint j = 0; j < ind; j++) nt += (*iter)[j];
        for (uint j = 0; j < ind; j++) nt += mutations[i][j];
        nt += nucleos_;
        // for (uint j = ind + 1; j < (*iter).nucleos.size(); j++) nt += (*iter)[j];
        for (uint j = ind + 1; j < mutations[i].nucleos.size(); j++) nt += mutations[i][j];
        // Update nucleos and size_modifier fields:
        // (*iter).nucleos = nt;
        mutations[i].nucleos = nt;
        // (*iter).size_modifier += size_;
        mutations[i].size_modifier += size_;
        // Adjust new positions and total scaffold size:
        // calc_positions(iter + 1, size_);
        calc_positions(i + 1, size_);
        ;
    /*
     If `new_pos_` is in the reference sequence following the Mutation, we add
     a new Mutation object:
     */
    } else {
        // uint p = iter - mutations.begin();
        ;
        // uint old_pos_ = ind + ((*iter).old_pos - (*iter).size_modifier);
        uint old_pos_ = ind + (mutations[i].old_pos - mutations[i].size_modifier);
        Mutation new_mut(old_pos_, new_pos_, nucleos_);
        // ++iter;
        ++i;
        mutations.insert(mutations.begin() + i, new_mut);
        // Adjust new positions and total scaffold size:
        // iter = mutations.begin() + p + 1;
        ;
        // calc_positions(iter);
        calc_positions(i);
    }
    return;
}





/*
 ------------------
 Add a substitution somewhere in the deque
 ------------------
 */
void VarSequence::add_substitution(const char& nucleo, const uint& new_pos_) {

    uint i = get_mut_ind_(new_pos_);

    // std::deque<Mutation>::iterator iter = get_mut_(new_pos_);
    // `mutations.end()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    // if (iter == mutations.end()) {
    if (i == mutations.size()) {
        std::string nucleos_(1, nucleo);
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nucleos_);
        mutations.push_front(new_mut);
        // if (!mutations.empty()) {
        //     // mutations[0] = new_mut;
        //     mutations.push_front(new_mut);
        //     // iter = mutations.begin();
        // } else {
        //     // ;
        //     mutations.push_back(new_mut);
        //     // iter = mutations.begin();
        //     // mutations.insert(iter, new_mut);
        // }
        // mutations.push_back(new_mut);
        // std::sort(mutations.begin(), mutations.end());
    } else {
        // uint ind = new_pos_ - (*iter).new_pos;
        uint ind = new_pos_ - mutations[i].new_pos;
        // If `new_pos_` is within the mutation sequence:
        // if (ind <= (*iter).size_modifier) {
        if (ind <= mutations[i].size_modifier) {
            // (*iter).nucleos[ind] = nucleo;
            mutations[i].nucleos[ind] = nucleo;
        // If `new_pos_` is in the reference sequence following the mutation:
        } else {
            // uint old_pos_ = ind + ((*iter).old_pos - (*iter).size_modifier);
            uint old_pos_ = ind + (mutations[i].old_pos - mutations[i].size_modifier);
            std::string nucleos_(1, nucleo);
            Mutation new_mut(old_pos_, new_pos_, nucleos_);
            // ++iter;
            ++i;
            // mutations.insert(iter, new_mut);
            mutations.insert(mutations.begin() + i, new_mut);
            // iter = mutations.begin();
        }
    }

    // std::string nucleos_(1, nucleo);
    // Mutation new_mut(new_pos_, new_pos_, nucleos_);
    // mutations.push_back(new_mut);

    return;
}



/*
 -------------------
 Internal function to "blowup" mutation(s) due to a deletion.
 By "blowup", I mean it removes substitutions and insertions if they're covered
 entirely by the deletion, and it merges any deletions that are contiguous.
 This function is designed to be used after `get_mut_`, and the output from that
 function should be the `iter` argument for this function.
 It is also designed for `calc_positions` to be used on `iter` afterward
 (bc this function alters `iter` to point to the position after the deletion),
 along with `size_mod` for the `modifier` argument
 (i.e., `calc_positions(iter, size_mod)`).
 Note that there's a check to ensure that this is never run when `mutations`
 is empty.
 */
void VarSequence::deletion_blowup_(
        std::deque<Mutation>::iterator& iter,
        uint& deletion_start, uint& deletion_end, sint& size_mod) {

    /*
     Difference between first and last objects to be removed.
     See the while loop below for more info.
     */
    sint n_iters = 0;

    /*
     `mutations.end()` is returned from the `get_mut_` function if
     `deletion_start` is before the first Mutation object.
     */
    if (iter == mutations.end()) {
        iter = mutations.begin();
    /*
     If it's a substitution and has a position < the deletion starting point,
     we can simply skip to the next one.
     If they're equal, we don't iterate.
     If it's > the deletion starting point, that should never happen if `get_mut_` is
     working properly, so we return an error.
     */
    } else if ((*iter).size_modifier == 0) {
        if ((*iter).new_pos < deletion_start) {
            ++iter;
        } else if ((*iter).new_pos == deletion_start) {
            ;
        } else {
            stop("Iterator problem in deletion_blowup_");
        }
    /*
     If the first Mutation is an insertion, we may have to merge it with
     this deletion, and this is done with `merge_del_ins_`.
     This function will iterate to the next Mutation (and should never cause
     iterator invalidation even after erasing an insertion).
     It also adjusts `size_mod` appropriately.
     (`n_iters` is used later and doesn't matter here)
     */
    } else if ((*iter).size_modifier > 0) {
        // merge_del_ins_(iter, deletion_start, deletion_end, size_mod, n_iters);
        size_mod = 0;
        return;
    /*
     If it's a deletion and next to the new deletion, we merge their information
     before removing the old mutation.
     (`remove_mutation_` automatically moves the iterator to the location
     after the original.)

     If it's not next to the new deletion, we just iterate to the next mutation.
     */
    } else {
        if ((*iter).new_pos == deletion_start) {
            deletion_start += (*iter).size_modifier;
            size_mod += (*iter).size_modifier;
            remove_mutation_(iter);
        } else ++iter;
    }

    /*
     If `iter` no longer overlaps this deletion or if the deletion is gone (bc it
     absorbed part/all of an insertion), return now
     */
    if ((*iter).new_pos > deletion_end || size_mod == 0) return;

    /*
     If there is overlap, then we delete a range of Mutation objects.
     `iter` will point to the object after the last to be erased.
     `iter - n_iters` will point to the first object to be erased.
     */
    n_iters = 0;

    while (iter != mutations.end()) {
        if ((*iter).new_pos > deletion_end) break;
        // For substitutions, do nothing before iterating
        if ((*iter).size_modifier == 0) {
            ++iter;
            ++n_iters;
        /*
         For insertions, run `merge_del_ins_` to make sure that...
             (1) any sequence not overlapping the deletion is kept
             (2) `size_mod` is adjusted properly
             (3) insertions that are entirely overlapped by the deletion are erased
             (4) `iter` is moved to the next Mutation
             (5) `n_iters` is increased if the insertion isn't deleted
         */
        } else if ((*iter).size_modifier > 0) {
            // merge_del_ins_(iter, deletion_start, deletion_end, size_mod, n_iters);
            // // as above, stop here if deletion is absorbed
            // if (size_mod == 0) return;
            size_mod = 0;
            return;
        // For deletions, merge them with the current one
        } else {
            deletion_end -= (*iter).size_modifier;
            deletion_end = std::min(deletion_end, seq_size - 1);
            size_mod = deletion_start - deletion_end - 1;
            ++iter;
            ++n_iters;
        }
    }
    if (n_iters > iter - mutations.begin()) stop("n_iters is counted wrong");
    auto iter_begin = iter - n_iters;
    // Remove all mutations in the specified range:
    remove_mutation_(iter_begin, iter);

    // `iter` now points to the position AFTER the erasing.
}





/*
 Inner function to merge an insertion and deletion.
 `iter` points to the focal insertion.
 Deletion start and end points are for the new, variant sequence.
 `size_mod` is the size_modifier field for the Mutation object that will be
 created for the deletion. I change this value by the number of "virtual" nucleotides
 removed during this operation. Virtual nucleotides are the extra ones stored
 in an insertion's Mutation object (i.e., the ones other than the reference sequence).
 `n_iters` is how many positions were moved during this function.
 It also moves the iterator to the next Mutation object.
 */
void VarSequence::merge_del_ins_(std::deque<Mutation>::iterator& insertion,
                                  uint& deletion_start, uint& deletion_end,
                                  sint& size_mod, sint& n_iters) {

    // The starting and ending positions of the focal insertion
    uint& insertion_start((*insertion).new_pos);
    uint insertion_end = insertion_start + (*insertion).size_modifier;

    /*
     If the deletion doesn't overlap, move to the next Mutation
     */
    if (deletion_start > insertion_end || deletion_end < insertion_start) {
        ++insertion;
        ++n_iters;
    /*
     Else if the entire insertion is covered by the deletion, adjust size_mod and
     remove the Mutation object for the insertion:
     */
    } else if (deletion_start <= insertion_start && deletion_end >= insertion_end) {
        size_mod += (*insertion).size_modifier; // making it less negative
        // `remove_mutation_` moves to the next Mutation, so no need to do it here.
        // See its code below for more info.
        remove_mutation_(insertion);
    /*
     Else if there is overlap, adjust the size_mod, remove that part of
     the inserted sequence, and adjust the insertion's size modifier:
     */
    } else {

        // index for first char to erase from `nucleos`
        sint tmp = deletion_start - insertion_start;
        tmp = std::max(0, tmp);
        uint erase_ind0 = static_cast<uint>(tmp);
        // index for last char NOT to erase from `nucleos`
        uint erase_ind1 = deletion_end - insertion_start + 1;
        erase_ind1 = std::min(
            erase_ind1,
            static_cast<uint>((*insertion).nucleos.size())
        );

        // Adjust the size modifier for the eventual Mutation object
        // for the deletion (making it less negative)
        size_mod += (erase_ind1 - erase_ind0);

        std::string& nts((*insertion).nucleos);
        nts.erase(nts.begin() + erase_ind0, nts.begin() + erase_ind1);
        // // clear memory:
        // std::string(nts.begin(), nts.end()).swap(nts);

        // Adjust the insertion's size modifier
        (*insertion).size_modifier = (*insertion).nucleos.size() - 1;

        /*
         If this deletion removes the first part of the insertion but doesn't reach
         the end of the insertion, we have to adjust this insertion's `new_pos` manually
         and not iterate.

         This is because this insertion's starting position is not affected by
         the positions within this insertion that were removed.
         (All subsequent mutations' starting positions are.)

         Also, iterating will cause this mutation to be deleted.
         */
        if (deletion_start <= insertion_start && deletion_end < insertion_end) {
            (*insertion).new_pos += (erase_ind1 - erase_ind0);
        } else {
            ++insertion;
            ++n_iters;
        }
    }

    return;
}






/*
 Inner function to remove Mutation and keep iterator from being invalidated.
 After this function, `iter` points to the next item or `mutations.end()`.
 If two iterators are provided, the range of mutations are removed and each
 iterator now points to directly outside the range that was removed.
 If the removal occurs at the beginning of the mutations deque, then
 `iter1 == mutations.begin() && iter2 == mutations.begin()`
 after this function is run.
 */
void VarSequence::remove_mutation_(std::deque<Mutation>::iterator& mutation) {
    if (mutation == mutations.end()) return;
    uint p = mutation - mutations.begin();
    // erase:
    mutations.erase(mutation);
    // // clear memory:
    // std::deque<Mutation>(mutations.begin(), mutations.end()).swap(mutations);
    // reset iterator:
    mutation = mutations.begin() + p;
    return;
}
void VarSequence::remove_mutation_(std::deque<Mutation>::iterator& mutation1,
                                    std::deque<Mutation>::iterator& mutation2) {

    uint p1 = mutation1 - mutations.begin();
    // erase range:
    mutations.erase(mutation1, mutation2);
    // // clear memory:
    // std::deque<Mutation>(mutations.begin(), mutations.end()).swap(mutations);
    // reset iterators:
    if (p1 > 0) {
        mutation1 = mutations.begin() + p1 - 1;
        mutation2 = mutation1 + 1;
    } else {
        mutation1 = mutations.begin();
        mutation2 = mutations.begin();
    }

    return;
}



/*
 ------------------
 Internal function for finding character of either mutation or reference
 given an index (in the "new", variant sequence) and a single Mutation object.
 This only works if you've already narrowed it down to the Mutation object
 that is directly previous to the index position.
 ------------------
 */
char VarSequence::get_char_(const uint& new_pos,
                             const std::deque<Mutation>::iterator& mut) const {
    const Mutation& m(*mut);
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




/*
 ------------------
 Inner function to return an iterator to the Mutation object nearest to
 (without being past) an input position on the "new", variant scaffold.
 If the input position is before the first Mutation object or if `mutations` is empty,
 this function returns `mutations.end()`.
 ------------------
 */
std::deque<Mutation>::iterator VarSequence::get_mut_(const uint& new_pos) const {

    std::deque<Mutation>::iterator iter;

    if (mutations.empty()) return mutations.end();

    if (new_pos >= seq_size) {
        stop(
            "new_pos should never be >= the scaffold size. "
            "Either re-calculate the scaffold size or closely examine new_pos."
            );
    }
    /*
     If new_pos is less than the position for the first mutation, we return
     mutations.end():
     */
    if (new_pos < mutations.front().new_pos) {
        return mutations.end();
    /*
     If the new_pos is greater than or equal to the position for the last
     mutation, we return the last Mutation:
     */
    } else if (new_pos >= mutations.back().new_pos) {
        return mutations.end() - 1;
    }
    // /*
    //  If not either of the above, then we will first try to guess the approximate
    //  position to minimize how many iterations we have to perform.
    //  */
    // uint ind_guess = mutations.size() * new_pos / seq_size;
    // iter = mutations.begin() + ind_guess;
    // /*
    //  If the current mutation comes LATER than the input new position (`new_pos`),
    //  iterate BACKWARD in the mutation deque until the current mutation is
    //  EARLIER THAN OR EQUAL TO the new position.
    //  */
    // if ((*iter).new_pos > new_pos) {
    //     while ((*iter).new_pos > new_pos) --iter;
    //     /*
    //      If the current mutation comes EARLIER than the input new position (`new_pos`),
    //      iterate FORWARD in the mutation deque until the current mutation is
    //      EARLIER THAN OR EQUAL TO the new position.
    //
    //      Then, if IT REACHES THE END OF THE MUTATION VECTOR OR the current mutation is
    //      NOT EQUAL TO the input new position, iterate back to the previous mutation.
    //      I do this latter part because I want the Mutation object I use to relate
    //      to the positions it belongs to and the reference sequences AFTER it.
    //
    //      If the current mutation is EQUAL TO the input new position, check to make sure
    //      that it's not a deletion immediately followed by another mutation.
    //      If it is, then go to the next mutation.
    //      If not, then leave `out` unchanged.
    //      (It's assumed that there aren't two deletions right after each other bc that's
    //      a silly thing to do.)
    //      */
    // } else if ((*iter).new_pos < new_pos) {
    //     while (iter != mutations.end() && (*iter).new_pos < new_pos) ++iter;
    //     if (iter == mutations.end() || (*iter).new_pos != new_pos) {
    //         --iter;
    //     } else if ((*iter).new_pos == new_pos && (*iter).nucleos == "") {
    //         auto next = iter + 1;
    //         if (next != mutations.end() && (*next).new_pos == (*iter).new_pos) {
    //             ++iter;
    //         }
    //     }
    //     /*
    //      If the current mutation is EQUAL TO the input new position, we just need to
    //      check to make sure that it's not a deletion immediately followed by another
    //      mutation, as above.
    //      */
    // } else if ((*iter).new_pos == new_pos && (*iter).nucleos == "") {
    //     auto next = iter + 1;
    //     if (next != mutations.end() && (*next).new_pos == (*iter).new_pos) {
    //         ++iter;
    //     }
    // }
    //
    // // This is included for troubleshooting. It should be removed afterward.
    // if ((*iter).new_pos > new_pos) {
    //     stop(
    //         "Function `get_mut_` is returning a Mutation object whose "
    //         "new position field is greater than the input new position."
    //     );
    // }

    // Move mutation to the proper spot
    while (iter != mutations.end()) {
        if ((*iter).new_pos > new_pos) break;
        ++iter;
    }
    --iter;


    return iter;
}


uint VarSequence::get_mut_ind_(const uint& new_pos) const {

    if (mutations.empty()) return mutations.size();

    if (new_pos >= seq_size) {
        stop(
            "new_pos should never be >= the scaffold size. "
            "Either re-calculate the scaffold size or closely examine new_pos."
        );
    }
    /*
    If new_pos is less than the position for the first mutation, we return
    mutations.end():
    */
    if (new_pos < mutations.front().new_pos) {
        return mutations.size();
    /*
    If the new_pos is greater than or equal to the position for the last
    mutation, we return the last Mutation:
    */
    } else if (new_pos >= mutations.back().new_pos) {
        return mutations.size() - 1;
    }

    uint i = 0;
    // Move mutation to the proper spot
    while (mutations[i].new_pos <= new_pos) ++i;
    // Because we want the last mutation that is <= new_pos
    --i;

    return i;
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


// For printing info on a set of variants
void VarSet::print() const noexcept {

    uint total_muts = 0;
    for (const VarGenome& vg : variants) {
        for (const VarSequence& vs : vg.var_genome) {
            total_muts += vs.mutations.size();
        }
    }

    int console_width = get_width();

    int n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 21) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Variants object >>" << std::endl;

    Rcout.imbue(std::locale(""));
    Rcout << "# Variants: " << size() << std::endl;
    Rcout << "# Mutations: " << total_muts << std::endl;
    Rcout << std::endl;

    n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 28) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Reference genome info: >>" << std::endl;
    reference.print();
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
void see_ref(SEXP ref_) {
    XPtr<RefGenome> ref(ref_);
    ref->print();
    return;
}

//[[Rcpp::export]]
std::string get_ref_seq(SEXP ref_, const uint& s) {
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
std::string get_ref_name(SEXP ref_, const uint& s) {
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
uint get_ref_seq_size(SEXP ref_, const uint& s) {
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
uint get_ref_n_scaff(SEXP ref_) {
    XPtr<RefGenome> ref_xptr(ref_);
    RefGenome& ref(*ref_xptr);
    uint out = ref.size();
    return out;
}


//' Make a VarSet object from a set of sequences and # variants
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_vars(const std::deque<std::string>& seqs, const uint& n_vars) {
    XPtr<VarSet> vset(new VarSet(seqs, n_vars), true);
    return vset;
}


//' Function to piece together the strings for all sequences in a VarGenome.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
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


//' Function to print info on a VarSet.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_vs(SEXP vs_) {
    XPtr<VarSet> vs(vs_);
    vs->print();
    return;
}


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
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
List see_mutations(SEXP vs_, const uint& v) {

    XPtr<VarSet> vs(vs_);
    VarGenome& vg((*vs)[v]);

    List out(vg.size());
    for (uint i = 0; i < vg.size(); i++) {
        const VarSequence& vs(vg.var_genome[i]);
        std::vector<sint> size_mod;
        std::vector<uint> old_pos;
        std::vector<uint> new_pos;
        std::vector<std::string> nucleos;
        for (auto iter = vs.mutations.begin(); iter != vs.mutations.end(); ++iter) {
            size_mod.push_back((*iter).size_modifier);
            old_pos.push_back((*iter).old_pos);
            new_pos.push_back((*iter).new_pos);
            nucleos.push_back((*iter).nucleos);
        }
        DataFrame mutations_i = DataFrame::create(
            _["size_mod"] = size_mod,
            _["old_pos"] = old_pos,
            _["new_pos"] = new_pos,
            _["nucleos"] = nucleos);
        out[i] = mutations_i;
    }
    return out;
}


//' See all scaffold sizes in a VarSet object.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
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

// //[[Rcpp::export]]
// void add_substitution(SEXP vs_, const uint& v,
//                       const uint& scaff,
//                       const char& nucleo,
//                       const uint& new_pos_) {
//     XPtr<VarSet> vset(vs_);
//     VarGenome& vg((*vset)[v]);
//     VarSequence& vs(vg[scaff]);
//     vs.add_substitution(nucleo, new_pos_);
//     return;
// }
// //[[Rcpp::export]]
// void add_insertion(SEXP vs_, const uint& v,
//                    const uint& scaff,
//                    const std::string& nucleos_,
//                    const uint& new_pos_) {
//     XPtr<VarSet> vset(vs_);
//     VarGenome& vg((*vset)[v]);
//     VarSequence& vs(vg[scaff]);
//     vs.add_insertion(nucleos_, new_pos_);
//     return;
// }
// //[[Rcpp::export]]
// void add_deletion(SEXP vs_, const uint& v,
//                   const uint& scaff,
//                   const uint& size_,
//                   const uint& new_pos_) {
//     XPtr<VarSet> vset(vs_);
//     VarGenome& vg((*vset)[v]);
//     VarSequence& vs(vg[scaff]);
//     vs.add_deletion(size_, new_pos_);
//     return;
// }


//' View the starting portion of a variant sequence.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::string see_start(SEXP vs_, const uint& v,
                      const uint& scaff,
                      const uint& size_) {
    XPtr<VarSet> vset(vs_);
    VarGenome& vg((*vset)[v]);
    VarSequence& vs(vg[scaff]);
    std::string out = vs.get_seq_start(size_);
    return out;
}


//' View a chunk of a variant sequence.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
std::string see_chunk(SEXP vs_, const uint& v,
                      const uint& scaff,
                      const uint& start,
                      const uint& chunk_size) {
    XPtr<VarSet> vset(vs_);
    VarGenome& vg((*vset)[v]);
    VarSequence& vs(vg[scaff]);
    std::string out;
    std::deque<Mutation>::iterator mut = vs.mutations.begin();

    vs.set_seq_chunk(out, start, chunk_size, mut);

    return out;
}





// void many_mutations(const std::deque<std::string>& seqs, const uint& n_vars,


//' Add many mutations (> 1,000) to a VarSet object from R.
//'
//' I made this fxn for testing bc running add_* (which goes from R to C++) many
//' times gives segfault when including deletions.
//' `min_muts` and `max_muts` give range of # mutations per variant sequence.
//'
//' Temporary function for testing.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void many_mutations(SEXP vs_,
                    const double& min_muts, const double& max_muts) {

    XPtr<VarSet> vs_xptr(vs_);
    VarSet& vset(*vs_xptr);

    // VarSet vset(seqs, n_vars);

    double prev_type;

    for (uint v = 0; v < vset.size(); v++) {
        for (uint s = 0; s < vset.reference.size(); s++) {
            VarSequence& vs(vset[v][s]);
            uint n_muts = static_cast<uint>(R::runif(min_muts, max_muts+1));
            uint m = 0;
            while (m < n_muts) {
                uint max_size = vs.seq_size;
                uint pos = static_cast<uint>(R::unif_rand() *
                    static_cast<double>(max_size));
                double rnd = R::unif_rand();
                if (rnd < 0.5) {
                    std::string str = cpp_rando_seq(1);
                    vs.add_substitution(str[0], pos);
                // } else if (rnd < 0.75) {
                } else {
                    uint size = static_cast<uint>(R::rexp(2.0) + 1.0);
                    if (size > 10) size = 10;
                    std::string str = cpp_rando_seq(size + 1);
                    vs.add_insertion(str, pos);
                }
                // else {
                //     uint size = static_cast<uint>(R::rexp(2.0) + 1.0);
                //     if (size > 10) size = 10;
                //     // vs.add_deletion(size, pos);
                // }
                prev_type = rnd;
                ++m;
            }
        }
    }

    return;
}
