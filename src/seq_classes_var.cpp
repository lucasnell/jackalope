
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


#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "util.h"  // get_width


using namespace Rcpp;



VarSequence& VarSequence::operator+=(const VarSequence& other) {

    // If either is empty, then this is easy:
    if (other.mutations.empty()) return *this;
    if (mutations.empty()) {
        mutations = other.mutations;
        seq_size = other.seq_size;
        return *this;
    }

    // Combine sequence sizes:
    sint64 diff = static_cast<sint64>(other.seq_size) -
        static_cast<sint64>(ref_seq->size());
    seq_size += diff;

    /*
    Now combine `mutations` deques.
    Process differently depending on whether `other` has its mutations before
    or after this one's.
    If they overlap, then throw an error.
    */
    // `other` has mutations before this one:
    bool other_is_before = other.mutations.back() < mutations.front();
    // `other` has mutations after this one:
    bool other_is_after = other.mutations.front() > mutations.back();
    if (other_is_before) {
        // Adjust current mutations' `new_pos` fields (using `diff` from above):
        auto mut_ = mutations.begin();
        for (; mut_ != mutations.end(); ++mut_) {
            (*mut_).new_pos += diff;
        }
        // Now add the new mutations:
        auto mut = other.mutations.rbegin();  // note the use of reverse iterator!
        for (; mut != other.mutations.rend(); ++mut) {
            // Add the new mutation to the front of `(*this).mutations`:
            mutations.push_front(*mut);
        }
    } else if (other_is_after) {
        // The amount to adjust the new mutation's `new_pos` fields:
        diff = static_cast<sint64>(seq_size) - static_cast<sint64>(ref_seq->size());
        auto mut = other.mutations.begin();
        for (; mut != other.mutations.end(); ++mut) {
            // Add the new mutation to the back of `(*this).mutations`:
            mutations.push_back(*mut);
            // Adjust the `new_pos` field:
            mutations.back().new_pos += diff;
        }
    } else {
        str_stop({"\nOverlapping VarSequence.mutations in +=. ",
                 "Note that when combining VarSequence objects, you must ",
                 "do it sequentially, either from the front or back."});
    }

    return *this;
}


/*
 ------------------
 Retrieve a nucleotide (char type) from the variant sequence
 based on the position in the new, variant sequence
 ------------------
 */
char VarSequence::get_nt(const uint64& new_pos) const {
    char out;
    /*
     Index to the Mutation object nearest to (without being past)
     an input position (see below for the `get_mut_` fxn):
     */
    uint64 mut_i = get_mut_(new_pos);
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
    uint64 mut_i = 0;

    std::string out(seq_size, 'x');
    uint64 pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations[mut_i].new_pos) {
        out[pos] = (*ref_seq)[pos];
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint64 next_mut_i = mut_i + 1;
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
                                const uint64& start,
                                const uint64& chunk_size,
                                uint64& mut_i) const {

    uint64 end = start + chunk_size - 1;

    if (start >= seq_size) {
        mut_i = mutations.size();
        chunk_str.clear();
        return;
    }
    // Making sure end doesn't go beyond the sequence bounds
    if (end >= seq_size) end = seq_size - 1;

    uint64 out_length = end - start + 1;

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
    // Clearing string if necessary (reserving memory should happen outside this method)
    if (chunk_str.size() > 0) chunk_str.clear();

    uint64 pos = start;
    uint64 next_mut_i = mut_i + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `start` is before the first mutation)
     */
    while (pos < mutations[mut_i].new_pos && pos <= end) {
        chunk_str += (*ref_seq)[pos];
        ++pos;
    }
    if (pos > end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `end`) but before the next mutation
     */
    while (next_mut_i < mutations.size()) {
        while (pos < mutations[next_mut_i].new_pos && pos <= end) {
            chunk_str += get_char_(pos, mut_i);
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
        chunk_str += get_char_(pos, mut_i);
        ++pos;
    }

    return;
}




/*
 Similar to above, but it's for use in the sequence simulator, and it fills info for
 the read.
 An important difference is that it doesn't always fill from the start of the
 read, since I don't want to muck with the barcodes.
 */

void VarSequence::fill_read(std::string& read,
                            const uint64& read_start,
                            const uint64& seq_start,
                            uint64 n_to_add) const {

    uint64 mut_i = 0;

    uint64 seq_end = seq_start + n_to_add - 1;
    // Making sure seq_end doesn't go beyond the sequence bounds
    if (seq_end >= seq_size) {
        seq_end = seq_size - 1;
        n_to_add = seq_size - seq_start;
    }

    // Make sure the read is long enough (this fxn should never shorten it):
    if (read.size() < n_to_add + read_start) read.resize(n_to_add + read_start, 'N');

    // No need to mess around with mutations if there aren't any
    if (mutations.empty()) {
        for (uint64 i = 0; i < n_to_add; i++) {
            read[(read_start + i)] = ref_seq->nucleos[(seq_start + i)];
        }
        return;
    }
    // Move mutation to the proper spot
    while (mut_i < mutations.size()) {
        if (seq_start < mutations[mut_i].new_pos) break;
        ++mut_i;
    }
    if (mut_i != 0) --mut_i;

    uint64 seq_pos = seq_start;
    uint64 read_pos = read_start;
    uint64 next_mut_i = mut_i + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `seq_start` is before the first mutation)
     */
    while (seq_pos < mutations[mut_i].new_pos && seq_pos <= seq_end) {
        read[read_pos] = (*ref_seq)[seq_pos];
        ++seq_pos;
        ++read_pos;
    }
    if (seq_pos > seq_end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `seq_end`) but before the next mutation
     */
    while (next_mut_i < mutations.size()) {
        while (seq_pos < mutations[next_mut_i].new_pos && seq_pos <= seq_end) {
            read[read_pos] = get_char_(seq_pos, mut_i);
            ++seq_pos;
            ++read_pos;
        }
        if (seq_pos > seq_end) return;
        ++mut_i;
        ++next_mut_i;
    }

    /*
     If we reach the last mutation, add nucleotides until `seq_end` (remember that above,
     I've made sure that `seq_end < seq_size`).
     */
    while (seq_pos <= seq_end) {
        read[read_pos] = get_char_(seq_pos, mut_i);
        ++seq_pos;
        ++read_pos;
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
char VarSequence::get_char_(const uint64& new_pos,
                            const uint64& mut_i) const {
    const Mutation& m(mutations[mut_i]);
    char out;
    uint64 ind = new_pos - m.new_pos;
    if (static_cast<sint64>(ind) > m.size_modifier) {
        ind += (m.old_pos - m.size_modifier);
        out = (*ref_seq)[ind];
    } else {
        out = m[ind];
    }
    return out;
}





/*
 ------------------
 Re-calculate new positions (and total sequence size)
 ------------------
 */

/*
 For ALL Mutation objects (this is only used when reading VCF files)
 */
void VarSequence::calc_positions() {

    if (mutations.size() == 0) return;

    uint64 mut_i = 0;

    sint64 modifier = mutations[mut_i].size_modifier;
    ++mut_i;

    // Updating individual Mutation objects
    for (; mut_i < mutations.size(); mut_i++) {
        mutations[mut_i].new_pos += modifier;
        modifier += mutations[mut_i].size_modifier;
    }
    // Updating full sequence size
    seq_size += modifier;

    return;
}
/*
 For all Mutation objects after a given Mutation object
 (this is for after you insert a NEW Mutation, where `mut_i` below points to
 that Mutation)
 */
void VarSequence::calc_positions(uint64 mut_i) {

    sint64 modifier = mutations[mut_i].size_modifier;
    ++mut_i;

    // Updating individual Mutation objects
    for (; mut_i < mutations.size(); mut_i++) {
        mutations[mut_i].new_pos += modifier;
    }
    // Updating full sequence size
    seq_size += modifier;

    return;
}
/*
 For all Mutation objects after AND INCLUDING a given Mutation object
 (this is for after you MERGE multiple Mutations, where `mut_i` below points to
 that merged Mutation and `modifier` refers to the net change in sequence size
 after the merge)
 */
void VarSequence::calc_positions(uint64 mut_i, const sint64& modifier) {
    // Updating individual Mutation objects
    for (; mut_i < mutations.size(); ++mut_i) {
        mutations[mut_i].new_pos += modifier;
    }
    // Updating full sequence size
    seq_size += modifier;

    return;
}












/*
 ------------------
 Add a deletion somewhere in the deque
 ------------------
 */
void VarSequence::add_deletion(const uint64& size_, const uint64& new_pos_) {

    if (size_ == 0 || new_pos_ >= seq_size) return;

    uint64 mut_i;

    // Renaming this for a more descriptive name and to allow it to change
    uint64 deletion_start = new_pos_;

    /*
     Last position this deletion refers to
     (`std::min` is to coerce this deletion to one that's possible):
     */
    uint64 deletion_end = std::min(deletion_start + size_ - 1, seq_size - 1);

    /*
     Position on the old reference sequence for this deletion.
     This will change if there are insertions or deletions before `new_pos_`.
     */
    uint64 old_pos_ = new_pos_;

    /*
     Size modifier of this deletion. This can change when deletion merges with an
     insertion or another deletion.
     */
    sint64 size_mod = deletion_start - deletion_end - 1;

    /*
     If `mutations` is empty, just add to the beginning and adjust sequence size
     */
    if (mutations.empty()) {
        Mutation new_mut(old_pos_, deletion_start, size_mod);
        mutations.push_front(new_mut);
        seq_size += size_mod;

        /*
         If the mutations deque isn't empty, we may need to edit some of the mutations
         after this deletion if they're affected by it.
         See `deletion_blowup_` below for more info.
         */
    } else {

        /*
         Sequence-size modifier to be used to edit subsequent mutations.
         This number does not change.
         */
        const sint64 subseq_modifier(size_mod);

        mut_i = get_mut_(deletion_start);

        /*
         This "blows up" subsequent mutations if they're destroyed/altered
         by this deletion.
         See `deletion_blowup_` below for more info.
         */
        deletion_blowup_(mut_i, deletion_start, deletion_end, size_mod);

        /*
         If `size_mod` is zero, this means that an insertion/insertions absorbed all
         of the deletion, so after adjusting sizes, our business is done here.
         */
        if (size_mod == 0) {
            calc_positions(mut_i, subseq_modifier);
            return;
        }

        /*
         If the deletion hasn't been absorbed, we need to calculate its
         position on the old (i.e., reference) sequence:
         */
        if (mut_i != 0) {
            --mut_i;
            old_pos_ = deletion_start - mutations[mut_i].new_pos +
            mutations[mut_i].old_pos - mutations[mut_i].size_modifier;
            ++mut_i;
        } else old_pos_ = deletion_start; // (`deletion_start` may have changed)


        // Adjust (1) positions of all mutations after and including `mut_i`, and
        //        (2) the sequence size
        calc_positions(mut_i, subseq_modifier);

        // Now create the Mutation and insert it.
        Mutation new_mut(old_pos_, deletion_start, size_mod);
        mutations.insert(mutations.begin() + mut_i, new_mut);
    }
    return;
}




/*
 ------------------
 Add an insertion somewhere in the deque
 ------------------
 */
void VarSequence::add_insertion(const std::string& nucleos_, const uint64& new_pos_) {

    uint64 mut_i = get_mut_(new_pos_);
    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nt = (*ref_seq)[new_pos_] + nucleos_;
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nt);
        mutations.push_front(new_mut);
        // Adjust new positions and total sequence size:
        calc_positions(static_cast<uint64>(0));
        return;
    }

    uint64 ind = new_pos_ - mutations[mut_i].new_pos;
    /*
     If `new_pos_` is within the Mutation sequence (which is never the case for
     deletions), then we adjust it as such:
     */
    if (static_cast<sint64>(ind) <= mutations[mut_i].size_modifier) {
        sint64 size_ = nucleos_.size();
        // string to store combined nucleotides
        std::string nt = "";
        for (uint64 j = 0; j <= ind; j++) nt += mutations[mut_i][j];
        nt += nucleos_;
        for (uint64 j = ind + 1; j < mutations[mut_i].nucleos.size(); j++) {
            nt += mutations[mut_i][j];
        }
        // Update nucleos and size_modifier fields:
        mutations[mut_i].nucleos = nt;
        mutations[mut_i].size_modifier += size_;
        // Adjust new positions and total sequence size:
        calc_positions(mut_i + 1, size_);
        /*
         If `new_pos_` is in the reference sequence following the Mutation, we add
         a new Mutation object:
         */
    } else {
        uint64 old_pos_ = ind + (mutations[mut_i].old_pos -
            mutations[mut_i].size_modifier);
        std::string nt = (*ref_seq)[old_pos_] + nucleos_;
        Mutation new_mut(old_pos_, new_pos_, nt);
        ++mut_i;
        mutations.insert(mutations.begin() + mut_i, new_mut);
        // Adjust new positions and total sequence size:
        calc_positions(mut_i);
    }
    return;
}





/*
 ------------------
 Add a substitution somewhere in the deque
 ------------------
 */
void VarSequence::add_substitution(const char& nucleo, const uint64& new_pos_) {

    uint64 mut_i = get_mut_(new_pos_);

    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nucleos_(1, nucleo);
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nucleos_);
        mutations.push_front(new_mut);
    } else {
        uint64 ind = new_pos_ - mutations[mut_i].new_pos;
        // If `new_pos_` is within the mutation sequence:
        if (static_cast<sint64>(ind) <= mutations[mut_i].size_modifier) {
            mutations[mut_i].nucleos[ind] = nucleo;
            // If `new_pos_` is in the reference sequence following the mutation:
        } else {
            uint64 old_pos_ = ind + (mutations[mut_i].old_pos -
                mutations[mut_i].size_modifier);
            std::string nucleos_(1, nucleo);
            Mutation new_mut(old_pos_, new_pos_, nucleos_);
            ++mut_i;
            mutations.insert(mutations.begin() + mut_i, new_mut);
        }
    }

    return;
}



/*
 -------------------
 Internal function to "blowup" mutation(s) due to a deletion.
 By "blowup", I mean it removes substitutions and insertions if they're covered
 entirely by the deletion, and it merges any deletions that are contiguous.
 This function is designed to be used after `get_mut_`, and the output from that
 function should be the `mut_i` argument for this function.
 It is also designed for `calc_positions` to be used on `mut_i` afterward
 (bc this function alters `mut_i` to point to the position after the deletion),
 along with `size_mod` for the `modifier` argument
 (i.e., `calc_positions(mut_i, size_mod)`).
 Note that there's a check to ensure that this is never run when `mutations`
 is empty.
 */
void VarSequence::deletion_blowup_(uint64& mut_i, uint64& deletion_start,
                                   uint64& deletion_end, sint64& size_mod) {

    /*
     ---------
     Taking care of the initial mutation pointed to:
     ---------
     */

    /*
     `mutations.size()` is returned from the `get_mut_` function if
     `deletion_start` is before the first Mutation object.
     */
    if (mut_i == mutations.size()) {
        mut_i = 0;
        /*
         If it's a substitution and has a position < the deletion starting point,
         we can simply skip to the next one.
         If they're equal, we don't iterate.
         If it's > the deletion starting point, that should never happen if `get_mut_` is
         working properly, so we return an error.
         */
    } else if (mutations[mut_i].size_modifier == 0) {
        if (mutations[mut_i].new_pos < deletion_start) {
            ++mut_i;
        } else if (mutations[mut_i].new_pos == deletion_start) {
            ;
        } else {
            stop("Index problem in deletion_blowup_");
        }
        /*
         If the first Mutation is an insertion, we may have to merge it with
         this deletion, and this is done with `merge_del_ins_`.
         This function will iterate to the next Mutation.
         It also adjusts `size_mod` appropriately.
         */
    } else if (mutations[mut_i].size_modifier > 0) {
        merge_del_ins_(mut_i, deletion_start, deletion_end, size_mod);
        /*
         If it's a deletion and next to the new deletion, we merge their information
         before removing the old mutation.
         (`remove_mutation_` automatically moves the index to the location
         after the original.)

         If it's not next to the new deletion, we just iterate to the next mutation.
         */
    } else {
        if (mutations[mut_i].new_pos == deletion_start) {
            size_mod += mutations[mut_i].size_modifier;
            remove_mutation_(mut_i);
        } else ++mut_i;
    }


    /*
     ---------
     Taking care of subsequent mutations:
     ---------
     */

    /*
     If `mut_i` no longer overlaps this deletion or if the deletion is gone (bc it
     absorbed part/all of an insertion), return now.
     (The first check is to prevent segfault when two deletions are the first
     to be added to a mutations deque.)
     */
    if (mut_i < mutations.size()) {
        if (mutations[mut_i].new_pos > deletion_end || size_mod == 0) return;
    }

    /*
     If there is overlap, then we delete a range of Mutation objects.
     `mut_i` will point to the object after the last to be erased.
     `range_begin` will point to the first object to be erased.
     */
    uint64 range_begin = mut_i;
    while (mut_i < mutations.size()) {
        if (mutations[mut_i].new_pos > deletion_end) break;
        // For substitutions, do nothing before iterating
        if (mutations[mut_i].size_modifier == 0) {
            ++mut_i;
            /*
             For insertions, run `merge_del_ins_` to make sure that...
             (1) any sequence not overlapping the deletion is kept
             (2) `size_mod` is adjusted properly
             (3) insertions that are entirely overlapped by the deletion are erased
             (4) `mut_i` is moved to the next Mutation
             */
        } else if (mutations[mut_i].size_modifier > 0) {
            merge_del_ins_(mut_i, deletion_start, deletion_end, size_mod);
            // as above, stop here if deletion is absorbed
            if (size_mod == 0) return;
            /*
             For deletions, merge them with the current one
             */
        } else {
            size_mod += mutations[mut_i].size_modifier;
            ++mut_i;
        }
    }

    // Remove all mutations in the specified range:
    remove_mutation_(range_begin, mut_i);

    // `mut_i` now points to the position AFTER the erasing.

    return;
}





/*
 Inner function to merge an insertion and deletion.
 `insert_i` points to the focal insertion.
 Deletion start and end points are for the new, variant sequence.
 `size_mod` is the size_modifier field for the Mutation object that will be
 created for the deletion. I change this value by the number of "virtual" nucleotides
 removed during this operation. Virtual nucleotides are the extra ones stored
 in an insertion's Mutation object (i.e., the ones other than the reference sequence).
 `n_iters` is how many positions were moved during this function.
 It also moves the index to the next Mutation object.
 */
void VarSequence::merge_del_ins_(uint64& insert_i,
                                 uint64& deletion_start, uint64& deletion_end,
                                 sint64& size_mod) {

    // The starting and ending positions of the focal insertion
    uint64& insertion_start(mutations[insert_i].new_pos);
    uint64 insertion_end = insertion_start + mutations[insert_i].size_modifier;

    /*
     If the deletion doesn't overlap, move to the next Mutation
     */
    if (deletion_start > insertion_end || deletion_end < insertion_start) {
        ++insert_i;
        /*
         Else if the entire insertion is covered by the deletion, adjust size_mod and
         remove the Mutation object for the insertion:
         */
    } else if (deletion_start <= insertion_start && deletion_end >= insertion_end) {
        size_mod += mutations[insert_i].size_modifier; // making it less negative
        /*
         Because we're deleting a mutation, `insert_i` refers to the next
         object without us doing anything here.
         */
        remove_mutation_(insert_i);
        /*
         Else if there is overlap, adjust the size_mod, remove that part of
         the inserted sequence, and adjust the insertion's size modifier:
         */
    } else {

        // index for first char to erase from `nucleos`
        sint64 tmp = deletion_start - insertion_start;
        if (tmp < 0L) tmp = 0L;
        uint64 erase_ind0 = static_cast<uint64>(tmp);
        // index for last char NOT to erase from `nucleos`
        uint64 erase_ind1 = deletion_end - insertion_start + 1;
        erase_ind1 = std::min(
            erase_ind1,
            static_cast<uint64>(mutations[insert_i].nucleos.size())
        );

        // Adjust the size modifier for the eventual Mutation object
        // for the deletion (making it less negative)
        size_mod += (erase_ind1 - erase_ind0);

        std::string& nts(mutations[insert_i].nucleos);
        nts.erase(nts.begin() + erase_ind0, nts.begin() + erase_ind1);
        // clear memory:
        clear_memory<std::string>(nts);

        // Adjust the insertion's size modifier
        mutations[insert_i].size_modifier = mutations[insert_i].nucleos.size() - 1;

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
            mutations[insert_i].new_pos += (erase_ind1 - erase_ind0);
        } else {
            ++insert_i;
        }
    }

    return;
}






/*
 Inner function to remove Mutation.
 After this function, `mut_i` points to the next item or `mutations.size()`.
 If two indices are provided, the range of mutations are removed and each
 index now points to directly outside the range that was removed.
 If the removal occurs at the beginning of the mutations deque, then
 `mut_i == 0 && mut_i2 == 0` after this function is run.
 */
void VarSequence::remove_mutation_(uint64& mut_i) {
    if (mut_i == mutations.size()) return;
    // erase:
    mutations.erase(mutations.begin() + mut_i);
    // clear memory:
    clear_memory<std::deque<Mutation>>(mutations);
    return;
}
void VarSequence::remove_mutation_(uint64& mut_i1, uint64& mut_i2) {

    // erase range:
    mutations.erase(mutations.begin() + mut_i1, mutations.begin() + mut_i2);
    // clear memory:
    clear_memory<std::deque<Mutation>>(mutations);
    // reset indices:
    if (mut_i1 > 0) {
        mut_i2 = mut_i1;
        mut_i1--;
    } else {
        mut_i1 = 0;
        mut_i2 = 0;
    }

    return;
}






/*
 ------------------
 Inner function to return an index to the Mutation object nearest to
 (without being past) an input position on the "new", variant sequence.
 If the input position is before the first Mutation object or if `mutations` is empty,
 this function returns `mutations.end()`.
 ------------------
 */

uint64 VarSequence::get_mut_(const uint64& new_pos) const {

    uint64 mut_i = 0;

    if (mutations.empty()) return mutations.size();

    if (new_pos >= seq_size) {
        str_stop({"new_pos should never be >= the sequence size. ",
                 "Either re-calculate the sequence size or closely examine new_pos."});

    }
    /*
     If new_pos is less than the position for the first mutation, we return
     mutations.size():
     */
    if (new_pos < mutations.front().new_pos) return mutations.size();

    /*
     If the new_pos is greater than or equal to the position for the last
     mutation, we return the last Mutation:
     */
    if (new_pos >= mutations.back().new_pos) return mutations.size() - 1;

    /*
     If not either of the above, then we will first try to guess the approximate
     position to minimize how many iterations we have to perform.
     */
    mut_i = static_cast<double>(mutations.size() * new_pos) /
        static_cast<double>(seq_size);
    /*
     If the current mutation is not past `new_pos`, iterate until it is.
     I'm intentionally going past the mutation to make sure we're not getting a deletion
     immediately followed by another mutation.
     (We don't need to check for `mut_i` getting to the last index
     (`mutations.size() - 1`) because we've already checked for that situation above.)
     */
    while (mutations[mut_i].new_pos <= new_pos) ++mut_i;
    /*
     Now move mutation to the proper spot: the last mutation that is <= `new_pos`.
     */
    while (mutations[mut_i].new_pos > new_pos) --mut_i;

    return mut_i;
}








void VarSet::print() const noexcept {

    uint64 total_muts = 0;
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

    Rcout << "# Variants: " << big_int_format<uint64>(size()) << std::endl;
    Rcout << "# Mutations: " << big_int_format<uint64>(total_muts) << std::endl;
    Rcout << std::endl;

    n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 28) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Reference genome info: >>" << std::endl;
    reference->print();
}
