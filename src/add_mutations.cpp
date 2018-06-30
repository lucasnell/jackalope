
/*
 ********************************************************

 Methods for adding mutations to sequence classes.
 Note that this file is for just adding the mutations, not for anything related to
 determining where they should go.

 ********************************************************
 */

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <deque>  // deque class
#include <algorithm>  // sort


#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes


using namespace Rcpp;






/*
 ------------------
 Re-calculate new positions (and total sequence size)
 ------------------
 */
/*
 For all Mutation objects after a given Mutation object
 (this is for after you insert a NEW Mutation, where `mut_i` below points to
 that Mutation)
 */
void VarSequence::calc_positions(uint32 mut_i) {

    sint32 modifier = mutations[mut_i].size_modifier;
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
void VarSequence::calc_positions(uint32 mut_i, const sint32& modifier) {
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
void VarSequence::add_deletion(const uint32& size_, const uint32& new_pos_) {

    uint32 mut_i;

    // Renaming this for a more descriptive name and to allow it to change
    uint32 deletion_start = new_pos_;

    /*
     Last position this deletion refers to
     (`std::min` is to coerce this deletion to one that's possible):
     */
    uint32 deletion_end = std::min(deletion_start + size_ - 1, seq_size - 1);

    /*
     Position on the old reference sequence for this deletion.
     This will change if there are insertions or deletions before `new_pos_`.
     */
    uint32 old_pos_ = new_pos_;

    /*
     Size modifier of this deletion. This can change when deletion merges with an
     insertion or another deletion.
     */
    sint32 size_mod = deletion_start - deletion_end - 1;

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
        const sint32 subseq_modifier(size_mod);

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
void VarSequence::add_insertion(const std::string& nucleos_, const uint32& new_pos_) {

    uint32 mut_i = get_mut_(new_pos_);
    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nt = ref_seq[new_pos_] + nucleos_;
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nt);
        mutations.push_front(new_mut);
        // Adjust new positions and total sequence size:
        calc_positions(static_cast<uint32>(0));
        return;
    }

    uint32 ind = new_pos_ - mutations[mut_i].new_pos;
    /*
     If `new_pos_` is within the Mutation sequence (which is never the case for
     deletions), then we adjust it as such:
     */
    if (static_cast<sint32>(ind) <= mutations[mut_i].size_modifier) {
        sint32 size_ = nucleos_.size();
        // string to store combined nucleotides
        std::string nt = "";
        for (uint32 j = 0; j <= ind; j++) nt += mutations[mut_i][j];
        nt += nucleos_;
        for (uint32 j = ind + 1; j < mutations[mut_i].nucleos.size(); j++) {
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
        uint32 old_pos_ = ind + (mutations[mut_i].old_pos -
            mutations[mut_i].size_modifier);
        std::string nt = ref_seq[old_pos_] + nucleos_;
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
void VarSequence::add_substitution(const char& nucleo, const uint32& new_pos_) {

    uint32 mut_i = get_mut_(new_pos_);

    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nucleos_(1, nucleo);
        // (below, notice that new position and old position are the same)
        Mutation new_mut(new_pos_, new_pos_, nucleos_);
        mutations.push_front(new_mut);
    } else {
        uint32 ind = new_pos_ - mutations[mut_i].new_pos;
        // If `new_pos_` is within the mutation sequence:
        if (static_cast<sint32>(ind) <= mutations[mut_i].size_modifier) {
            mutations[mut_i].nucleos[ind] = nucleo;
        // If `new_pos_` is in the reference sequence following the mutation:
        } else {
            uint32 old_pos_ = ind + (mutations[mut_i].old_pos -
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
void VarSequence::deletion_blowup_(uint32& mut_i, uint32& deletion_start,
                                   uint32& deletion_end, sint32& size_mod) {

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
     absorbed part/all of an insertion), return now
     */
    if (mutations[mut_i].new_pos > deletion_end || size_mod == 0) {
        return;
    }

    /*
     If there is overlap, then we delete a range of Mutation objects.
     `mut_i` will point to the object after the last to be erased.
     `range_begin` will point to the first object to be erased.
     */
    uint32 range_begin = mut_i;
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
void VarSequence::merge_del_ins_(uint32& insert_i,
                                 uint32& deletion_start, uint32& deletion_end,
                                 sint32& size_mod) {

    // The starting and ending positions of the focal insertion
    uint32& insertion_start(mutations[insert_i].new_pos);
    uint32 insertion_end = insertion_start + mutations[insert_i].size_modifier;

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
        sint32 tmp = deletion_start - insertion_start;
        if (tmp < 0L) tmp = 0L;
        uint32 erase_ind0 = static_cast<uint32>(tmp);
        // index for last char NOT to erase from `nucleos`
        uint32 erase_ind1 = deletion_end - insertion_start + 1;
        erase_ind1 = std::min(
            erase_ind1,
            static_cast<uint32>(mutations[insert_i].nucleos.size())
        );

        // Adjust the size modifier for the eventual Mutation object
        // for the deletion (making it less negative)
        size_mod += (erase_ind1 - erase_ind0);

        std::string& nts(mutations[insert_i].nucleos);
        nts.erase(nts.begin() + erase_ind0, nts.begin() + erase_ind1);
        // clear memory:
        std::string(nts.begin(), nts.end()).swap(nts);

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
void VarSequence::remove_mutation_(uint32& mut_i) {
    if (mut_i == mutations.size()) return;
    // erase:
    mutations.erase(mutations.begin() + mut_i);
    // clear memory:
    std::deque<Mutation>(mutations.begin(), mutations.end()).swap(mutations);
    return;
}
void VarSequence::remove_mutation_(uint32& mut_i1, uint32& mut_i2) {

    // erase range:
    mutations.erase(mutations.begin() + mut_i1, mutations.begin() + mut_i2);
    // clear memory:
    std::deque<Mutation>(mutations.begin(), mutations.end()).swap(mutations);
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

uint32 VarSequence::get_mut_(const uint32& new_pos) const {

    uint32 mut_i = 0;

    if (mutations.empty()) return mutations.size();

    if (new_pos >= seq_size) {
        stop(
            "new_pos should never be >= the sequence size. "
            "Either re-calculate the sequence size or closely examine new_pos."
            );
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



