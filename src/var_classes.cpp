
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
#include "hap_classes.h"  // Hap* classes
#include "util.h"  // get_width


using namespace Rcpp;







// `start` is inclusive
// this `HapChrom` must be empty after `mut_i`
// return `sint64` is the size modifier for mutations added
sint64 HapChrom::add_to_back(const HapChrom& other, const uint64& mut_i) {

    if (other.mutations.size() <= mut_i) return 0;

    if (!mutations.empty() &&
        mutations.old_pos.back() >= other.mutations.old_pos[mut_i]) {
        str_stop({"\nOverlapping HapChrom.mutations in HapChrom::add_to_back. ",
                 "Note that when combining HapChrom objects using `add_to_back`, you ",
                 "must do it sequentially, from the back ONLY."});
    }

    sint64 new_size_mod = 0;
    sint64 old_size_mod = static_cast<sint64>(chrom_size) -
        static_cast<sint64>(ref_chrom->size());

    for (uint64 i = mut_i; i < other.mutations.size(); i++) {
        mutations.push_back(other.mutations.old_pos[i],
                            other.mutations.new_pos[i],
                            other.mutations.nucleos[i]);
        mutations.new_pos.back() = mutations.old_pos.back() +
            old_size_mod + new_size_mod;
        new_size_mod += other.size_modifier(i);
    }

    chrom_size += new_size_mod;

    return new_size_mod;

}






/*
 ------------------
 Retrieve all nucleotides (i.e., the full chromosome; std::string type) from
 the haplotype chromosome
 ------------------
 */

std::string HapChrom::get_chrom_full() const {

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
 Set an input string object to any chunk of a chromosome from the haplotype chromosome.
 Before anything, this function moves `mut` to the location right before this chunk's
 starting position. I keep this index around so I don't have to iterate through
 the entire mutation deque multiple times.
 If end position is beyond the size of the chromosome, it changes `chunk_str` to the
 chromosome from the start to the chromosome end.
 If start position is beyond the size of the chromosome, it sets `mut` to
 `mutations.end()` and clears `chunk_str`.
 ------------------
 */
void HapChrom::set_chrom_chunk(std::string& chunk_str,
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

void HapChrom::fill_read(std::string& read,
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
void HapChrom::add_deletion(const uint64& size_, const uint64& new_pos_) {

    if (size_ == 0 || new_pos_ >= chrom_size) return;

    // Renaming this for a more descriptive name and to allow it to change
    uint64 deletion_start = new_pos_;
    uint64 deletion_end = std::min(deletion_start + size_ - 1, chrom_size - 1);
    const sint64 size_mod(deletion_start - deletion_end - 1);

    // If `mutations` is empty, just add to the beginning and adjust chromosome size
    if (mutations.empty()) {

        mutations.push_front(new_pos_, deletion_start, nullptr);
        chrom_size += size_mod;
        return;
    }
    /*
     If the first mutation is after the deletion,
     just add to the beginning, adjust the `new_pos`s, and adjust chromosome size
     */
    if (mutations.new_pos.front() > deletion_end) {

        /*
         In the rare case where the first mutation is a deletion that's right after
         the new deletion, we don't need to add an extra mutation but we do need to
         adjust the first mutation's `old_pos`.
         */
        bool del_after = mutations.new_pos.front() == (deletion_end + 1) &&
            size_modifier(0) < 0;

        for (uint64& np : mutations.new_pos) np += size_mod; // <-- gets done either way

        if (del_after) {
            mutations.old_pos.front() += size_mod;
        } else {
            mutations.push_front(new_pos_, deletion_start, nullptr);
        }
        chrom_size += size_mod;
        return;
    }

    /*
     This is for when the first mutation starts after the deletion but will be at
     least partially removed by it.
     This is a weird situation that needs to be addressed explicitly.
     */
    bool first_overlap = mutations.new_pos.front() > deletion_start &&
        mutations.new_pos.front() <= deletion_end;

    /*
     (Not using `get_mut_` below bc we want the first mutation that's == `deletion_start`
      or, if that doesn't exist, the last one that's < `deletion_start`.)
     */
    uint64 mut_i;
    if (mutations.new_pos.back() < deletion_start) {
        mut_i = mutations.size() - 1;
    } else {
        mut_i = 0;
        // Iterate until the first mutation that's >= `deletion_start`
        while (mutations.new_pos[mut_i] < deletion_start) ++mut_i;
        // Go back one if it's > `deletion_start` (But not if `mut_i` is zero!)
        if (mutations.new_pos[mut_i] > deletion_start && mut_i > 0) --mut_i;
    }

    // Getting old position info before changing any mutation info
    uint64 old_pos_ = deletion_old_pos_(deletion_start, deletion_end, mut_i);


    /*
     This (1) returns which mutations get removed because of this deletion,
     (2) adjusts size_mod for both insertions and (contiguous) deletions, and
     (3) edits insertions that are partially deleted.
     */
    std::vector<uint64> rm_inds;
    sint64 size_mod_remaining = size_mod; // to keep track of how much deletion remains
    for (uint64 i = mut_i; i < mutations.size(); i++) {
        deletion_one_mut_(i, deletion_start, deletion_end, size_mod,
                          size_mod_remaining, rm_inds);
    }

    chrom_size += size_mod;

    // Remove all mutations that have been deleted / merged:
    if (rm_inds.size() == 1) {
        mutations.erase(rm_inds.front());
    } else if (rm_inds.size() > 1) {
        mutations.erase(rm_inds.front(), rm_inds.back() + 1);
    } else {
        // Because in situations above, `rm_inds.front()` points to position after removal
        rm_inds = {mut_i};
        // Because in the first overlap scenario, we want to add it to the beginning:
        if (!first_overlap) rm_inds.front()++;
    }

    /*
     If `size_mod` is >= zero, this means that insertion(s) absorbed all
     of the deletion, so our business is done here.
     */
    if (size_mod_remaining >= 0) return;


    // Otherwise insert mutation info:
    mutations.insert(rm_inds.front(), old_pos_, deletion_start, nullptr);

    return;
}






/*
 ------------------
 Add an insertion somewhere in the deque
 ------------------
 */
void HapChrom::add_insertion(const std::string& nucleos_, const uint64& new_pos_) {

    sint64 size_mod = nucleos_.size();

    uint64 mut_i = get_mut_(new_pos_);
    // `mutations.size()` is returned above if `new_pos_` is before the
    // first mutation  or if `mutations` is empty
    if (mut_i == mutations.size()) {
        std::string nts = (*ref_chrom)[new_pos_] + nucleos_;
        // (below, notice that new position and old position are the same)
        mutations.push_front(new_pos_, new_pos_, nts.c_str());
        // Adjust new positions and total chromosome size:
        calc_positions(1, size_mod);
        return;
    }

    uint64 ind = new_pos_ - mutations.new_pos[mut_i];
    /*
     If `new_pos_` is within the Mutation chromosome (which is never the case for
     deletions), then we adjust it as such:
     */
    if (static_cast<sint64>(ind) <= size_modifier(mut_i)) {
        // string to store combined nucleotides
        std::string nts = "";
        for (uint64 j = 0; j <= ind; j++) nts += mutations.nucleos[mut_i][j];
        nts += nucleos_;
        uint64 nucleos_size = std::strlen(mutations.nucleos[mut_i]);
        for (uint64 j = ind + 1; j < nucleos_size; j++) {
            nts += mutations.nucleos[mut_i][j];
        }
        // Update nucleos field:
        delete [] mutations.nucleos[mut_i]; // delete old char array
        mutations.nucleos[mut_i] = new char[nts.size() + 1];
        std::copy(nts.begin(), nts.end(), mutations.nucleos[mut_i]);
        mutations.nucleos[mut_i][nts.size()] = '\0';
        // Adjust new positions and total chromosome size:
        calc_positions(mut_i + 1, size_mod);
        /*
         If `new_pos_` is in the reference chromosome following the Mutation, we add
         a new Mutation object:
         */
    } else {
        uint64 old_pos_ = ind + (mutations.old_pos[mut_i] -
            size_modifier(mut_i));
        std::string nts = (*ref_chrom)[old_pos_] + nucleos_;
        ++mut_i;
        mutations.insert(mut_i, old_pos_, new_pos_, nts.c_str());
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
void HapChrom::add_substitution(const char& nucleo, const uint64& new_pos_) {

    uint64 mut_i = get_mut_(new_pos_);

    // `mutations.size()` is returned above if `new_pos_` is before the
    // first Mutation object or if `mutations` is empty
    if (mut_i == mutations.size()) {
        // (below, notice that new position and old position are the same)
        mutations.push_front(new_pos_, new_pos_, nucleo);
    } else {
        uint64 ind = new_pos_ - mutations.new_pos[mut_i];
        // If `new_pos_` is within the mutation chromosome:
        if (static_cast<sint64>(ind) <= size_modifier(mut_i)) {
            /*
             If this new mutation reverts a substitution back to reference state,
             delete the Mutation object from the `mutations` field.
             Otherwise, adjust the mutation's sequence.
             */
            if ((size_modifier(mut_i) == 0) &&
                (ref_chrom->nucleos[mutations.old_pos[mut_i]] == nucleo)) {
                mutations.erase(mut_i);
            } else mutations.nucleos[mut_i][ind] = nucleo;
            // If `new_pos_` is in the reference chromosome following the mutation:
        } else {
            uint64 old_pos_ = ind + (mutations.old_pos[mut_i] -
                size_modifier(mut_i));
            ++mut_i;
            mutations.insert(mut_i, old_pos_, new_pos_, nucleo);
        }
    }

    return;
}






/*
 -------------------
 Inner function to get old position for deletion.
 -------------------
 */
uint64 HapChrom::deletion_old_pos_(const uint64& deletion_start,
                                    const uint64& deletion_end,
                                    const uint64& mut_i) const {

    if (mutations.new_pos[mut_i] == deletion_start) {
        return mutations.old_pos[mut_i];
    }
    /*
     This is for when the first mutation starts after the deletion but will be at
     least partially removed by it.
     */
    if (mutations.new_pos[mut_i] > deletion_start) {
        return deletion_start;
    }

    sint64 sm = size_modifier(mut_i);
    // (below can overflow if sm < 0, but that's fine bc it won't be used in that case.)
    uint64 mut_end = mutations.new_pos[mut_i] + sm;

    // This works only for subs and deletions, plus for insertions that aren't overlapping
    if (sm <= 0 || mut_end < deletion_start) {
        uint64 old_pos = deletion_start - mutations.new_pos[mut_i] +
            mutations.old_pos[mut_i] - sm;
        return old_pos;
    }

    /*
     This is true if the deletion...
     (1) takes out the first chunk of the insertion but some inserted sequence remains
     (2) takes out the whole insertion.
     For (1), this value won't be used bc no extra mutation will be added.
     For (2), this is the right value to use for the new mutation.
     */
    if (deletion_start == mutations.new_pos[mut_i]) {
        return mutations.old_pos[mut_i];
    }


    /*
     This is true if the deletion...
     (1) takes out a middle / ending chunk of the insertion but some inserted
         sequence remains
     (2) takes out an ending chunk of the insertion but some deletion remains
     For (1), this value won't be used bc no extra mutation will be added.
     For (2), this is the right value to use for the new mutation.
     */
    return mutations.old_pos[mut_i] + 1;

}



/*
 -------------------
 Inner function to adjust a single mutation for a deletion.

 This (1) adds indices for mutations that get removed because of this deletion,
 (2) merges any deletions that are contiguous,
 (3) edits insertions that are partially deleted, and
 (4) adjusts `new_size_mod` to inform whether it's been absorbed by insertion(s).
 -------------------
 */
void HapChrom::deletion_one_mut_(const uint64& mut_i,
                                 const uint64& deletion_start,
                                 const uint64& deletion_end,
                                 const sint64& full_size_mod,
                                 sint64& new_size_mod,
                                 std::vector<uint64>& rm_inds) {

    uint64& mut_pos(mutations.new_pos[mut_i]);

    // If it's after (and not next to) the deletion, then adjust the new_pos and finish
    if (mut_pos > (deletion_end + 1)) {
        mut_pos += full_size_mod;
        return;
    }

    sint64 sm = size_modifier(mut_i);


    /*
     Substitutions
     */
    if (sm == 0) {
        // If it's immediately after the deletion, adjust new_pos and finish
        if (mut_pos > deletion_end) {
            mut_pos += full_size_mod;
            return;
        }
        // If it's before the deletion, do nothing
        if (mut_pos < deletion_start) return;
        /*
         If neither of the above is true, there must be overlap, so we mark this
         mutation for removal.
         */
        rm_inds.push_back(mut_i);
        return;
    }

    /*
     Insertions
     */
    if (sm > 0) {

        // If it's immediately after the deletion, adjust new_pos and finish
        if (mut_pos > deletion_end) {
            mut_pos += full_size_mod;
            return;
        }

        uint64 mut_end = mut_pos + sm;

        if (mut_end < deletion_start) return;


        /*
         If the entire insertion is covered by the deletion, adjust new_size_mod and
         mark this mutation for removal.
         */
        if (deletion_start <= mut_pos && deletion_end >= mut_end) {
            new_size_mod += sm; // making it less negative
            rm_inds.push_back(mut_i);
            return;
        }

        /*
         Else if there is overlap, adjust the new_size_mod and remove that part of
         the inserted chromosome:
         */

        // index for first char to erase from `nucleos`
        sint64 tmp = deletion_start - mut_pos;
        if (tmp < 0L) tmp = 0L;
        uint64 erase_ind0 = static_cast<uint64>(tmp);
        // index for last char NOT to erase from `nucleos`
        uint64 erase_ind1 = deletion_end - mut_pos + 1;
        erase_ind1 = std::min(
            erase_ind1,
            static_cast<uint64>(std::strlen(mutations.nucleos[mut_i]))
        );
        new_size_mod += (erase_ind1 - erase_ind0);

        /*
         Re-size nucleotides for this mutation.
         I'm doing this by making a std::string, resizing, then assinging it back
         to the `char*` object in `nucleos`
         */
        std::string nts(mutations.nucleos[mut_i]);
        nts.erase(nts.begin() + erase_ind0, nts.begin() + erase_ind1);
        // delete and re-assign:
        delete [] mutations.nucleos[mut_i];
        mutations.nucleos[mut_i] = new char[nts.size() + 1];
        std::strcpy(mutations.nucleos[mut_i], nts.c_str());


        /*
         Only adjust this insertion's `new_pos` if the deletion...
         (1) removes the first part of the insertion
         (2) does NOT reach the end of the insertion
         (3) has a starting position before this insertion
         */
        if (deletion_start < mut_pos && deletion_end < mut_end) {
            // mut_pos -= (erase_ind1 + erase_ind0);
            mut_pos += (erase_ind1 - erase_ind0);
            mut_pos += full_size_mod;
        }

        return;

    }

    /*
     Deletions
     */

    // If it's before the deletion, do nothing
    if (mut_pos < deletion_start) return;

    /*
     In any other situation, we merge the deletions' information and mark the
     old one for removal.
     */
    new_size_mod += sm; // making it more negative
    rm_inds.push_back(mut_i);

    return;
}











/*
 ------------------
 Inner function to return an index to the Mutation object nearest to
 (without being past) an input position on the "new", haplotype chromosome.
 If the input position is before the first Mutation object or if `mutations` is empty,
 this function returns `mutations.end()`.
 ------------------
 */

uint64 HapChrom::get_mut_(const uint64& new_pos) const {

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








void HapSet::print() const noexcept {

    uint64 total_muts = 0;
    for (const HapGenome& vg : haplotypes) {
        for (const HapChrom& vc : vg.chromosomes) {
            total_muts += vc.mutations.size();
        }
    }

    int console_width = get_width();

    int n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 21) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< haplotypes object >>" << std::endl;

    Rcout << "# Haplotypes: " << big_int_format<uint64>(size()) << std::endl;
    Rcout << "# Mutations: " << big_int_format<uint64>(total_muts) << std::endl;
    Rcout << std::endl;

    n_spaces = static_cast<int>(
        std::ceil(static_cast<double>(console_width - 28) / 2)
    );

    for (int i = 0; i < n_spaces; i++) Rcout << ' ';
    Rcout << "<< Reference genome info: >>" << std::endl;
    reference->print();
}
