
/*
 ********************************************************

 Basic methods for chromosome classes: retrieving info and initializing.

 ********************************************************
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <cstring>  // C strings, including std::strcpy
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque


#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "var_classes.h"  // Var* classes
#include "util.h"  // get_width


using namespace Rcpp;


// `start` is inclusive
// this `VarChrom` must be empty after `mut_i`
// return `sint64` is the size modifier for mutations added
sint64 VarChrom::add_to_back(const VarChrom& other, const uint64& mut_i) {

    if (other.mutations.size() <= mut_i) return 0;

    if (!mutations.empty() &&
        mutations.old_pos.back() >= other.mutations.old_pos[mut_i]) {
        str_stop({"\nOverlapping VarChrom.mutations in VarChrom::add_to_back. ",
                 "Note that when combining VarChrom objects using `add_to_back`, you ",
                 "must do it sequentially, from the back ONLY."});
    }

    sint64 new_size_mod = 0;
    sint64 old_size_mod = static_cast<sint64>(chrom_size) -
        static_cast<sint64>(ref_chrom->size());

    for (uint64 i = mut_i; i < other.mutations.size(); i++) {
        mutations.push_back(other.mutations.size_modifier[i],
                            other.mutations.old_pos[i],
                            other.mutations.new_pos[i],
                            other.mutations.nucleos[i]);
        mutations.new_pos.back() = mutations.old_pos.back() +
            old_size_mod + new_size_mod;
        new_size_mod += other.mutations.size_modifier[i];
    }

    chrom_size += new_size_mod;

    return new_size_mod;

}






/*
 ------------------
 Retrieve all nucleotides (i.e., the full chromosome; std::string type) from
 the variant chromosome
 ------------------
 */

std::string VarChrom::get_chrom_full() const {

    if (mutations.empty()) return ref_chrom->nucleos;

    // Index to the first Mutation object
    uint64 mut_i = 0;

    std::string out;
    out.reserve(chrom_size);
    uint64 pos = 0;

    // Picking up any nucleotides before the first mutation
    while (pos < mutations.new_pos[mut_i]) {
        out.push_back((*ref_chrom)[pos]);
        ++pos;
    }

    // Now, for each subsequent mutation except the last, add all nucleotides
    // at or after its position but before the next one
    uint64 next_mut_i = mut_i + 1;
    while (next_mut_i < mutations.size()) {
        while (pos < mutations.new_pos[next_mut_i]) {
            out.push_back(get_char_(pos, mut_i));
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos < chrom_size) {
        out.push_back(get_char_(pos, mut_i));
        ++pos;
    }

    return out;
}





/*
 ------------------
 Set an input string object to any chunk of a chromosome from the variant chromosome.
 Before anything, this function moves `mut` to the location right before this chunk's
 starting position. I keep this index around so I don't have to iterate through
 the entire mutation deque multiple times.
 If end position is beyond the size of the chromosome, it changes `chunk_str` to the
 chromosome from the start to the chromosome end.
 If start position is beyond the size of the chromosome, it sets `mut` to
 `mutations.end()` and clears `chunk_str`.
 ------------------
 */
void VarChrom::set_chrom_chunk(std::string& chunk_str,
                                const uint64& start,
                                const uint64& chunk_size,
                                uint64& mut_i) const {

    uint64 end = start + chunk_size - 1;

    if (start >= chrom_size) {
        mut_i = mutations.size();
        chunk_str.clear();
        return;
    }
    // Making sure end doesn't go beyond the chromosome bounds
    if (end >= chrom_size) end = chrom_size - 1;

    uint64 out_length = end - start + 1;

    // No need to mess around with mutations if there aren't any
    if (mutations.empty()) {
        chunk_str = ref_chrom->nucleos.substr(start, out_length);
        return;
    }
    // Move mutation to the proper spot
    while (mut_i < mutations.size()) {
        if (start < mutations.new_pos[mut_i]) break;
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
    while (pos < mutations.new_pos[mut_i] && pos <= end) {
        chunk_str += (*ref_chrom)[pos];
        ++pos;
    }
    if (pos > end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `end`) but before the next mutation
     */
    while (next_mut_i < mutations.size()) {
        while (pos < mutations.new_pos[next_mut_i] && pos <= end) {
            chunk_str += get_char_(pos, mut_i);
            ++pos;
        }
        if (pos > end) return;
        ++mut_i;
        ++next_mut_i;
    }

    /*
     If we reach the last mutation, add nucleotides until `end` (remember that above,
     I've made sure that `end < chrom_size`).
     */
    while (pos <= end) {
        chunk_str += get_char_(pos, mut_i);
        ++pos;
    }

    return;
}




/*
 Similar to above, but it's for use in the chromosome simulator, and it fills info for
 the read.
 An important difference is that it doesn't always fill from the start of the
 read, since I don't want to muck with the barcodes.
 */

void VarChrom::fill_read(std::string& read,
                            const uint64& read_start,
                            const uint64& chrom_start,
                            uint64 n_to_add) const {

    uint64 mut_i = 0;

    uint64 chrom_end = chrom_start + n_to_add - 1;
    // Making sure chrom_end doesn't go beyond the chromosome bounds
    if (chrom_end >= chrom_size) {
        chrom_end = chrom_size - 1;
        n_to_add = chrom_size - chrom_start;
    }

    // Make sure the read is long enough (this fxn should never shorten it):
    if (read.size() < n_to_add + read_start) read.resize(n_to_add + read_start, 'N');

    // No need to mess around with mutations if there aren't any
    if (mutations.empty()) {
        for (uint64 i = 0; i < n_to_add; i++) {
            read[(read_start + i)] = ref_chrom->nucleos[(chrom_start + i)];
        }
        return;
    }
    // Move mutation to the proper spot
    while (mut_i < mutations.size()) {
        if (chrom_start < mutations.new_pos[mut_i]) break;
        ++mut_i;
    }
    if (mut_i != 0) --mut_i;

    uint64 chrom_pos = chrom_start;
    uint64 read_pos = read_start;
    uint64 next_mut_i = mut_i + 1;

    /*
     Picking up any nucleotides before the focal mutation (this should only happen when
     `mut == mutations.begin()` and `chrom_start` is before the first mutation)
     */
    while (chrom_pos < mutations.new_pos[mut_i] && chrom_pos <= chrom_end) {
        read[read_pos] = (*ref_chrom)[chrom_pos];
        ++chrom_pos;
        ++read_pos;
    }
    if (chrom_pos > chrom_end) return;

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position (and `chrom_end`) but before the next mutation
     */
    while (next_mut_i < mutations.size()) {
        while (chrom_pos < mutations.new_pos[next_mut_i] && chrom_pos <= chrom_end) {
            read[read_pos] = get_char_(chrom_pos, mut_i);
            ++chrom_pos;
            ++read_pos;
        }
        if (chrom_pos > chrom_end) return;
        ++mut_i;
        ++next_mut_i;
    }

    /*
     If we reach the last mutation, add nucleotides until `chrom_end`
     (remember that above, I've made sure that `chrom_end < chrom_size`).
     */
    while (chrom_pos <= chrom_end) {
        read[read_pos] = get_char_(chrom_pos, mut_i);
        ++chrom_pos;
        ++read_pos;
    }

    return;
}






/*
 ------------------
 Add a deletion somewhere in the deque
 ------------------
 */
void VarChrom::add_deletion(const uint64& size_, const uint64& new_pos_) {

    if (size_ == 0 || new_pos_ >= chrom_size) return;

    uint64 mut_i;

    // Renaming this for a more descriptive name and to allow it to change
    uint64 deletion_start = new_pos_;

    /*
     Last position this deletion refers to
     (`std::min` is to coerce this deletion to one that's possible):
     */
    uint64 deletion_end = std::min(deletion_start + size_ - 1, chrom_size - 1);

    /*
     Position on the old reference chromosome for this deletion.
     This will change if there are insertions or deletions before `new_pos_`.
     */
    uint64 old_pos_ = new_pos_;

    /*
     Size modifier of this deletion. This can change when deletion merges with an
     insertion or another deletion.
     */
    sint64 size_mod = deletion_start - deletion_end - 1;

    /*
     If `mutations` is empty, just add to the beginning and adjust chromosome size
     */
    if (mutations.empty()) {

        mutations.push_front(size_mod, old_pos_, deletion_start, nullptr);
        chrom_size += size_mod;

        /*
         If the mutations deque isn't empty, we may need to edit some of the mutations
         after this deletion if they're affected by it.
         See `deletion_blowup_` below for more info.
         */
    } else {

        /*
         Chromosome-size modifier to be used to edit subsequent mutations.
         This number does not change.
         */
        const sint64 subchrom_modifier(size_mod);

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
            calc_positions(mut_i, subchrom_modifier);
            return;
        }

        /*
         If the deletion hasn't been absorbed, we need to calculate its
         position on the old (i.e., reference) chromosome:
         */
        if (mut_i != 0) {
            --mut_i;
            old_pos_ = deletion_start - mutations.new_pos[mut_i] +
                mutations.old_pos[mut_i] - mutations.size_modifier[mut_i];
            ++mut_i;
        } else old_pos_ = deletion_start; // (`deletion_start` may have changed)


        // Adjust (1) positions of all mutations after and including `mut_i`, and
        //        (2) the chromosome size
        calc_positions(mut_i, subchrom_modifier);

        // Now insert mutation info:
        mutations.insert(mut_i, size_mod, old_pos_, deletion_start, nullptr);
    }
    return;
}




/*
 ------------------
 Add an insertion somewhere in the deque
 ------------------
 */
void VarChrom::add_insertion(const std::string& nucleos_, const uint64& new_pos_) {

    sint64 size_mod = nucleos_.size();

    uint64 mut_i = get_mut_(new_pos_);
    // `mutations.size()` is returned above if `new_pos_` is before the
    // first mutation  or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nts = (*ref_chrom)[new_pos_] + nucleos_;
        // (below, notice that new position and old position are the same)
        mutations.push_front(size_mod, new_pos_, new_pos_, nts.c_str());
        // Adjust new positions and total chromosome size:
        calc_positions(1, size_mod);
        return;
    }

    uint64 ind = new_pos_ - mutations.new_pos[mut_i];
    /*
     If `new_pos_` is within the Mutation chromosome (which is never the case for
     deletions), then we adjust it as such:
     */
    if (static_cast<sint64>(ind) <= mutations.size_modifier[mut_i]) {
        // string to store combined nucleotides
        std::string nts = "";
        for (uint64 j = 0; j <= ind; j++) nts += mutations.nucleos[mut_i][j];
        nts += nucleos_;
        uint64 nucleos_size = std::strlen(mutations.nucleos[mut_i]);
        for (uint64 j = ind + 1; j < nucleos_size; j++) {
            nts += mutations.nucleos[mut_i][j];
        }
        // Update nucleos and size_modifier fields:
        delete [] mutations.nucleos[mut_i]; // delete old char array
        mutations.nucleos[mut_i] = new char[nts.size() + 1];
        std::copy(nts.begin(), nts.end(), mutations.nucleos[mut_i]);
        mutations.nucleos[mut_i][nts.size()] = '\0';
        mutations.size_modifier[mut_i] += size_mod;
        // Adjust new positions and total chromosome size:
        calc_positions(mut_i + 1, size_mod);
        /*
         If `new_pos_` is in the reference chromosome following the Mutation, we add
         a new Mutation object:
         */
    } else {
        uint64 old_pos_ = ind + (mutations.old_pos[mut_i] -
            mutations.size_modifier[mut_i]);
        std::string nts = (*ref_chrom)[old_pos_] + nucleos_;
        ++mut_i;
        mutations.insert(mut_i, size_mod, old_pos_, new_pos_, nts.c_str());
        // Adjust new positions and total chromosome size:
        calc_positions(mut_i + 1, size_mod);
    }
    return;
}





/*
 ------------------
 Add a substitution somewhere in the deque
 ------------------
 */
void VarChrom::add_substitution(const char& nucleo, const uint64& new_pos_) {

    uint64 mut_i = get_mut_(new_pos_);

    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        // (below, notice that new position and old position are the same)
        mutations.push_front(0, new_pos_, new_pos_, nucleo);
    } else {
        uint64 ind = new_pos_ - mutations.new_pos[mut_i];
        // If `new_pos_` is within the mutation chromosome:
        if (static_cast<sint64>(ind) <= mutations.size_modifier[mut_i]) {
            /*
             If this new mutation reverts a substitution back to reference state,
             delete the Mutation object from the `mutations` field.
             Otherwise, adjust the mutation's sequence.
             */
            if ((mutations.size_modifier[mut_i] == 0) &&
                (ref_chrom->nucleos[mutations.old_pos[mut_i]] == nucleo)) {
                mutations.erase(mut_i);
            } else mutations.nucleos[mut_i][ind] = nucleo;
            // If `new_pos_` is in the reference chromosome following the mutation:
        } else {
            uint64 old_pos_ = ind + (mutations.old_pos[mut_i] -
                mutations.size_modifier[mut_i]);
            ++mut_i;
            mutations.insert(mut_i, 0, old_pos_, new_pos_, nucleo);
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
void VarChrom::deletion_blowup_(uint64& mut_i,
                                uint64& deletion_start,
                                uint64& deletion_end,
                                sint64& size_mod) {

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
    } else if (mutations.size_modifier[mut_i] == 0) {
        if (mutations.new_pos[mut_i] < deletion_start) {
            ++mut_i;
        } else if (mutations.new_pos[mut_i] == deletion_start) {
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
    } else if (mutations.size_modifier[mut_i] > 0) {
        merge_del_ins_(mut_i, deletion_start, deletion_end, size_mod);
        /*
         If it's a deletion and next to the new deletion, we merge their information
         before removing the old mutation.
         (`remove_mutation_` automatically moves the index to the location
         after the original.)

         If it's not next to the new deletion, we just iterate to the next mutation.
         */
    } else {
        if (mutations.new_pos[mut_i] == deletion_start) {
            size_mod += mutations.size_modifier[mut_i];
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
        if (mutations.new_pos[mut_i] > deletion_end || size_mod == 0) return;
    }

    /*
     If there is overlap, then we delete a range of Mutation objects.
     `mut_i` will point to the object after the last to be erased.
     `range_begin` will point to the first object to be erased.
     */
    uint64 range_begin = mut_i;
    while (mut_i < mutations.size()) {
        if (mutations.new_pos[mut_i] > deletion_end) break;
        // For substitutions, do nothing before iterating
        if (mutations.size_modifier[mut_i] == 0) {
            ++mut_i;
            /*
             For insertions, run `merge_del_ins_` to make sure that...
             (1) any chromosome not overlapping the deletion is kept
             (2) `size_mod` is adjusted properly
             (3) insertions that are entirely overlapped by the deletion are erased
             (4) `mut_i` is moved to the next Mutation
             */
        } else if (mutations.size_modifier[mut_i] > 0) {
            merge_del_ins_(mut_i, deletion_start, deletion_end, size_mod);
            // as above, stop here if deletion is absorbed
            if (size_mod == 0) return;
            /*
             For deletions, merge them with the current one
             */
        } else {
            size_mod += mutations.size_modifier[mut_i];
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
 Deletion start and end points are for the new, variant chromosome.
 `size_mod` is the size_modifier field for the Mutation object that will be
 created for the deletion. I change this value by the number of "virtual" nucleotides
 removed during this operation. Virtual nucleotides are the extra ones stored
 in an insertion's Mutation object (i.e., the ones other than the reference chromosome).
 `n_iters` is how many positions were moved during this function.
 It also moves the index to the next Mutation object.
 */
void VarChrom::merge_del_ins_(uint64& insert_i,
                              uint64& deletion_start,
                              uint64& deletion_end,
                              sint64& size_mod) {

    // The starting and ending positions of the focal insertion
    uint64& insertion_start(mutations.new_pos[insert_i]);
    uint64 insertion_end = insertion_start + mutations.size_modifier[insert_i];

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
        size_mod += mutations.size_modifier[insert_i]; // making it less negative
        /*
         Because we're deleting a mutation, `insert_i` refers to the next
         object without us doing anything here.
         */
        remove_mutation_(insert_i);
        /*
         Else if there is overlap, adjust the size_mod, remove that part of
         the inserted chromosome, and adjust the insertion's size modifier:
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
            static_cast<uint64>(std::strlen(mutations.nucleos[insert_i]))
        );

        // Adjust the size modifier for the eventual Mutation object
        // for the deletion (making it less negative)
        size_mod += (erase_ind1 - erase_ind0);

        /*
         Re-size nucleotides for this mutation.
         I'm doing this by making a std::string, resizing, then assinging it back
         to the `char*` object in `nucleos`
         */
        std::string nts(mutations.nucleos[insert_i]);
        nts.erase(nts.begin() + erase_ind0, nts.begin() + erase_ind1);
        // delete and re-assign:
        delete [] mutations.nucleos[insert_i];
        mutations.nucleos[insert_i] = new char[nts.size() + 1];
        std::strcpy(mutations.nucleos[insert_i], nts.c_str());

        // Adjust the insertion's size modifier
        mutations.size_modifier[insert_i] = static_cast<sint64>(nts.size()) - 1;

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
            mutations.new_pos[insert_i] += (erase_ind1 - erase_ind0);
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
void VarChrom::remove_mutation_(uint64& mut_i) {
    if (mut_i == mutations.size()) return;
    // erase:
    mutations.erase(mut_i);
    return;
}
void VarChrom::remove_mutation_(uint64& mut_i1, uint64& mut_i2) {

    // erase range:
    mutations.erase(mut_i1, mut_i2);

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
 (without being past) an input position on the "new", variant chromosome.
 If the input position is before the first Mutation object or if `mutations` is empty,
 this function returns `mutations.end()`.
 ------------------
 */

uint64 VarChrom::get_mut_(const uint64& new_pos) const {

    uint64 mut_i = 0;

    if (mutations.empty()) return mutations.size();

    if (new_pos >= chrom_size) {
        str_stop({"new_pos should never be >= the chromosome size. ",
                 "Either re-calculate the chromosome size or closely examine new_pos."});

    }
    /*
     If new_pos is less than the position for the first mutation, we return
     mutations.size():
     */
    if (new_pos < mutations.new_pos.front()) return mutations.size();

    /*
     If the new_pos is greater than or equal to the position for the last
     mutation, we return the last Mutation:
     */
    if (new_pos >= mutations.new_pos.back()) return mutations.size() - 1;

    /*
     If not either of the above, then we will first try to guess the approximate
     position to minimize how many iterations we have to perform.
     */
    mut_i = static_cast<double>(mutations.size() * new_pos) /
        static_cast<double>(chrom_size);
    /*
     If the current mutation is not past `new_pos`, iterate until it is.
     I'm intentionally going past the mutation to make sure we're not getting a deletion
     immediately followed by another mutation.
     (We don't need to check for `mut_i` getting to the last index
     (`mutations.size() - 1`) because we've already checked for that situation above.)
     */
    while (mutations.new_pos[mut_i] <= new_pos) ++mut_i;
    /*
     Now move mutation to the proper spot: the last mutation that is <= `new_pos`.
     */
    while (mutations.new_pos[mut_i] > new_pos) --mut_i;

    return mut_i;
}








void VarSet::print() const noexcept {

    uint64 total_muts = 0;
    for (const VarGenome& vg : variants) {
        for (const VarChrom& vc : vg.var_genome) {
            total_muts += vc.mutations.size();
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
