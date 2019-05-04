
/*
 This defines classes for sampling mutation locations with weights based on
 the mutation rate at each position.
 */

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "util.h"  // str_stop
#include "mutator_location.h"



using namespace Rcpp;




/*
 Adjust for a deletion.
 If the deletion totally overlaps it, this function changes this regions's rate
 to -1, and makes both start and end 0;
 it lastly returns true.
 Otherwise, it returns false.
*/
bool GammaRegion::deletion_adjust(const uint32& del_start,
                                  const uint32& del_end,
                                  const uint32& del_size) {

    if (rate == -1) return false;

    // No overlap and deletion starts after it
    if (del_start > end) return false;
    // No overlap and deletion starts before it
    if (del_end < start) {
        start -= del_size;
        end -= del_size;
        return false;
    }
    /*
     ----- Total overlap  -----
    */
    if ((del_start <= start) && (del_end >= end)) {
        start = 0;
        end = 0;
        rate = -1;
        return true;
    }
    /*
     Deletion is totally inside this region but doesn't entirely overlap it
     or touch any of the end points
     */
    if ((del_start > start) && (del_end < end)) {
        end -= del_size;
        return false;
    }
    // Partial overlap at the start
    if ((del_end >= start) && (del_start <= start)) {
        start = del_start;
        end -= del_size;
        return false;
    }
    // Partial overlap at the end
    if ((del_start <= end) && (del_end >= end)) {
        end = del_start - 1;
    }
    return false;
}





/*
 Process one row in Gamma matrix, adding GammaRegion(s) associated with it to the
 `regions` field.
 It also splits up each GammaRegion so that it only ever refers to a region of
 size `gamma_size`, or as close to this as possible.
 */
inline void LocationSampler::one_gamma_row(const arma::mat& gamma_mat,
                                           const uint32& i,
                                           uint32& mut_i,
                                           std::vector<uint32>& sizes) {

    GammaRegion reg;
    reg.gamma = gamma_mat(i,1); // this will be same even if it gets split

    uint32 start, size, n_regs;
    if (i > 0) {
        /*
         Even though `gamma_mat.col(0)` uses 1-based indices, I'm not subtracting 1
         below because I'm looking at previous row:
         */
        start = static_cast<uint32>(gamma_mat(i-1,0));
    } else start = 0;

    /*
     Because `size = end - start + 1` and bc `end = gamma_mat(i,0) - 1`
     (`gamma_mat.col(0)` uses 1-based indices),
     we don't need to add or subtract 1.
     */
    size = static_cast<uint32>(gamma_mat(i,0)) - start;

    n_regs = std::round(static_cast<double>(size) / static_cast<double>(gamma_size));

    // Vector of # bases per GammaRegion
    if (n_regs > 1) {
        sizes = split_int(size, n_regs);
    } else {
        sizes = std::vector<uint32>(1, size);
    }

    // String to store sequence info:
    std::string seq;
    seq.reserve(sizes[0] + 1);

    for (uint32 j = 0; j < sizes.size(); j++) {

        // Set bounds:
        reg.start = start;
        reg.end = start + sizes[j] - 1;

        // Calculate rate:
        reg.rate = 0;
        // (in `set_seq_chunk` below, `seq` gets cleared before it's filled)
        var_seq->set_seq_chunk(seq, reg.start, sizes[j], mut_i);
        for (const char& c : seq) reg.rate += nt_rates[c];
        reg.rate *= reg.gamma;

        // Set LocationSampler info:
        total_rate += reg.rate;
        regions.push_back(reg);

        // Update start for next iteration:
        start += sizes[j];
    }

    return;
}





void LocationSampler::construct_gammas(arma::mat gamma_mat) {

    if (gamma_size < 1) stop("Gamma size cannot be < 1.");

    total_rate = 0;

    if (!var_seq) stop("Cannot do construct_gammas method when var_seq isn't set.");
    // Sort from first to last region
    arma::uvec sort_inds = arma::sort_index(gamma_mat.col(0));
    gamma_mat = gamma_mat.rows(sort_inds);
    if (gamma_mat(gamma_mat.n_rows - 1, 0) != var_seq->size()) {
        stop("gamma_mat doesn't end with size of variant sequence");
    }
    // Now clear and fill in the regions vector
    regions.clear();
    regions.reserve((var_seq->size() / gamma_size) * 1.5);

    uint32 mut_i = 0;
    std::vector<uint32> sizes;
    sizes.reserve(1000); // arbitrarily chosen

    for (uint32 i = 0; i < gamma_mat.n_rows; i++) {
        one_gamma_row(gamma_mat, i, mut_i, sizes);
    }

    clear_memory<std::vector<GammaRegion>>(regions);

    end_rate = total_rate;

    return;
}



void LocationSampler::update_gamma_regions(const sint32& size_change,
                                           const uint32& pos) {

    /*
    -----------
    Substitutions
    -----------
    */
    if (size_change == 0) return;


    /*
    -----------
    InDels
    -----------
    */

    uint32 idx = get_gamma_idx(pos);


    /*
    Insertions
    */
    if (size_change > 0) {
        regions[idx].end += size_change;
        idx++;
        // update all following ranges:
        while (idx < regions.size()) {
            if (regions[idx].rate >= 0) {
                regions[idx].end += size_change;
                regions[idx].start += size_change;
            }
            idx++;
        }

        return;
    }

    /*
     Deletions
     */
    const uint32& del_start(pos);
    uint32 del_size = std::abs(size_change);
    uint32 del_end = pos + del_size - 1;

    // Iterate through and adjust all regions including and following the deletion:
    // std::vector<uint32> erase_inds;
    while (idx < regions.size()) {
        if (regions[idx].deletion_adjust(del_start, del_end, del_size)) {
            // erase_inds.push_back(idx);
        }
        idx++;
    }

    // /*
    // If any regions need erasing, their indices will be stored in erase_inds.
    // They will be consecutive indices, so we only need to access the front and back.
    // */
    // if (erase_inds.size() == 1) {
    //     regions.erase(regions.begin() + erase_inds.front());
    // } else if (erase_inds.size() > 1) {
    //     regions.erase(regions.begin() + erase_inds.front(),
    //                   regions.begin() + erase_inds.back() + 1);
    // }

    return;
}





// Used to check if gamma region needs to be iterated to the next one.
inline void LocationSampler::check_gamma(const uint32& pos,
                        uint32& gamma_end,
                        uint32& gam_i,
                        double& gamma) const {

    if (pos > gamma_end) {
        gam_i++;
        // Skip over deleted region(s) if necessary:
        while (gam_i < regions.size() && regions[gam_i].rate < 0) gam_i++;
        // Now set info:
        gamma = regions[gam_i].gamma;
        gamma_end = regions[gam_i].end;
    }
    return;
}




// To return the overall rate for an entire sequence or part of it:

// Inner method that does most of the work for `calc_rate` below

double LocationSampler::calc_rate__(uint32 start, uint32 end) const {

    double out = 0;

    if (var_seq->size() == 0) return out;

    /*
     If there are no mutations or if `end` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     (I'm using separate statements to avoid calling `front()` on an empty deque.)
     */
    bool use_mutations = true;
    if (var_seq->mutations.empty()) {
        use_mutations = false;
    } else if (var_seq->mutations.front().new_pos > end) {
        use_mutations = false;
    }
    if (!use_mutations) {

        uint32 i = start, gam_i = get_gamma_idx(start);

        while (i <= end) {
            double gamma = regions[gam_i].gamma;
            double tmp = 0;
            while (i <= regions[gam_i].end && i <= end) {
                tmp += nt_rates[var_seq->ref_seq->nucleos[i]];
                i++;
            }
            out += (tmp * gamma);
            gam_i++;
            // If the next one's been deleted, then keep going:
            while (gam_i < regions.size() && regions[gam_i].rate < 0) gam_i++;
        }

        return out;
    }


    // Index to the first Mutation object not past `start` position:
    uint32 mut_i = var_seq->get_mut_(start);
    // Index to the corresponding gamma region:
    uint32 gam_i = get_gamma_idx(start);

    double gamma = regions[gam_i].gamma;
    uint32 gamma_end = regions[gam_i].end;

    // Current position
    uint32 pos = start;

    /*
     If `start` is before the first mutation (resulting in
     `mut_i == var_seq->mutations.size()`),
     we must pick up any nucleotides before the first mutation.
     */
    if (mut_i == var_seq->mutations.size()) {
        mut_i = 0;
        for (; pos < var_seq->mutations[mut_i].new_pos; pos++) {
            check_gamma(pos, gamma_end, gam_i, gamma);
            out += (nt_rates[(*(var_seq->ref_seq))[pos]] * gamma);
        }
        check_gamma(pos, gamma_end, gam_i, gamma);
    }


    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position but before the next one.
     I'm adding `pos <= end` inside all while-statement checks to make sure
     it doesn't keep going after we've reached `end`.
     */
    uint32 next_mut_i = mut_i + 1;
    while (pos <= end && next_mut_i < var_seq->mutations.size()) {
        while (pos <= end && pos < var_seq->mutations[next_mut_i].new_pos) {
            char c = var_seq->get_char_(pos, mut_i);
            out += nt_rates[c] * gamma;
            ++pos;
            check_gamma(pos, gamma_end, gam_i, gamma);
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos <= end &&pos < var_seq->seq_size) {
        char c = var_seq->get_char_(pos, mut_i);
        out += nt_rates[c] * gamma;
        ++pos;
        check_gamma(pos, gamma_end, gam_i, gamma);
    }

    return out;
}



double LocationSampler::calc_rate() const {

    if (var_seq->size() == 0) return 0.0;

    uint32 start = 0;
    uint32 end = var_seq->size() - 1;

    double out = calc_rate__(start, end);

    return out;
}

double LocationSampler::calc_rate(const uint32& start, const uint32& end) const {

    if (var_seq->size() == 0) return 0.0;

    double out = calc_rate__(start, end);

    return out;
}








double LocationSampler::deletion_rate_change(const uint32& del_size,
                                             const uint32& start) {

    uint32 end = start + del_size - 1;
    if (end >= var_seq->size()) end = var_seq->size() - 1;

    std::string seq;
    seq.reserve(del_size);
    uint32 mut_i;
    safe_get_mut(start, mut_i);

    var_seq->set_seq_chunk(seq, start, end - start + 1, mut_i);

    double r, out = 0;
    uint32 seq_i = 0;
    uint32 idx = get_gamma_idx(start);

    while (seq_i < seq.size()) {
        GammaRegion& reg(regions[idx]);
        // Skip over deleted region if necessary:
        if (reg.rate < 0) {
            idx++;
            continue;
        }
        while (((seq_i + start) <= reg.end) && (seq_i < seq.size())) {
            r = reg.gamma * nt_rates[seq[seq_i]];
            out -= r;
            reg.rate -= r;
            total_rate -= r;
            end_rate -= r;
            seq_i++;
        }
        idx++;
    }

    return out;
}







/*
 This more-efficiently retrieves the index to the Mutation object inside mutations
 field, much like `VarSequence::get_mut_`.
 The difference is that it returns 0 in the situations where the`get_mut_` returns
 mutations.size(). This works more safely for methods here.
 */
inline void LocationSampler::safe_get_mut(const uint32& pos, uint32& mut_i) const {

    mut_i = 0;

    const std::deque<Mutation>& mutations(var_seq->mutations);

    /*
     If new_pos is less than the position for the first mutation, we return
     keep it at zero.
    */
    if (mutations.empty() || (pos < mutations.front().new_pos)) return;

    /*
     If the new_pos is greater than or equal to the position for the last
     mutation, we return the last Mutation:
     */
    if (pos >= mutations.back().new_pos) {
        mut_i = mutations.size() - 1;
        return;
    }

    /*
     If not either of the above, then we will first try to guess the approximate
     position to minimize how many iterations we have to perform.
     */
    mut_i = static_cast<double>(mutations.size() * pos) /
        static_cast<double>(var_seq->size());
    /*
     If the current mutation is not past `pos`, iterate until it is.
     I'm intentionally going past the mutation to make sure we're not getting a deletion
     immediately followed by another mutation.
     (We don't need to check for `mut_i` getting to the last index
     (`mutations.size() - 1`) because we've already checked for that situation above.)
     */
    while (mutations[mut_i].new_pos <= pos) ++mut_i;
    /*
     Now move mutation to the proper spot: the last mutation that is <= `pos`.
     */
    while (mutations[mut_i].new_pos > pos) --mut_i;

    return;

}





/*
 For a given position within a gamma region, find the rate associated with it.
 It just returns the rate (inclusive) from `reg.start` to `end`.
 This function should never be fed a deleted region!
 */
inline long double LocationSampler::partial_gamma_rate___(
        const uint32& end,
        const GammaRegion& reg) const {

    long double out = 0;

    if (end > reg.end) stop("end > reg.end");
    if (end < reg.start) stop("end < reg.start");
    if (reg.rate < 0) stop("partial_gamma_rate___ should not be run on a deleted region");
    if (end == reg.end || reg.rate == 0) return reg.rate;

    /*
    If there are no mutations or if `end` is before the first mutation,
    then we don't need to use the `mutations` field at all.
    */
    if (var_seq->mutations.empty() || (end < var_seq->mutations.front().new_pos)) {

        for (uint32 i = reg.start; i <= end; i++) {
            out += nt_rates[var_seq->ref_seq->nucleos[i]];
        }
        out *= reg.gamma;
        return out;

    }

    // Current position
    uint32 pos = reg.start;
    // Index to the first Mutation object not past `pos`:
    uint32 mut_i = var_seq->get_mut_(pos);

    /*
     If `pos` is before the first mutation (resulting in
     `mut_i == var_seq->mutations.size()`),
     we must pick up any nucleotides before the first mutation.
     */
    if (mut_i == var_seq->mutations.size()) {
        mut_i = 0;
        for (; pos < var_seq->mutations[mut_i].new_pos; pos++) {
            out += nt_rates[var_seq->ref_seq->nucleos[pos]];
        }
    }

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position but before the next one.
     I'm adding `pos <= end` inside all while-statement checks to make sure
     it doesn't keep going after we've reached `end`.
     */
    uint32 next_mut_i = mut_i + 1;
    while (pos <= end && next_mut_i < var_seq->mutations.size()) {
        while (pos <= end && pos < var_seq->mutations[next_mut_i].new_pos) {
            char c = var_seq->get_char_(pos, mut_i);
            out += nt_rates[c];
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos <= end && pos < var_seq->seq_size) {
        char c = var_seq->get_char_(pos, mut_i);
        out += nt_rates[c];
        ++pos;
    }

    out *= reg.gamma;

    return out;

}



/*
 Update starting and ending positions and rates
 */
void LocationSampler::update_start_end(const uint32& start, const uint32& end) {

    // If they don't need updated, then don't do it:
    if (start_end_set && (start == start_pos) && (end == end_pos)) return;

    if ((start >= var_seq->size()) || (end >= var_seq->size())) {
        stop("start or end in update_start_end is >= variant sequence size");
    }
    if (start > end) stop("start > end in update_start_end");

    start_pos = start;
    end_pos = end;

    // If they point to full sequence starts and ends, then this is easy:
    if ((start == 0) && (end == var_seq->size())) {
        start_rate = 0;
        end_rate = total_rate;
        return;
    }

    long double cum_wt = 0;
    uint32 gamm_i = 0;
    // Skip over deleted region(s) if necessary:
    while (regions[gamm_i].rate < 0) gamm_i++;
    long double part_rate = 0;

    uint32 mut_i = 0;

    /*
     Rate for starting position:
     */
    if (start_pos == 0) {
        start_rate = 0;
    } else {
        /*
         Below, I'm using `start_pos - 1` for things because to sample for positions
         from `start` to `end` inclusively, I want the starting rate to be for the
         position before `start`
         */
        safe_get_mut(start_pos - 1, mut_i);
        cum_wt = regions[gamm_i].rate;
        // Find the GammaRegion for the new start:
        while ((gamm_i < (regions.size() - 1)) &&
               (regions[gamm_i].end < (start_pos - 1))) {
            gamm_i++;
            if (regions[gamm_i].rate > 0) cum_wt += regions[gamm_i].rate;
        }
        if (regions[gamm_i].rate < 0) stop("update should not end on deleted region");
        cum_wt -= regions[gamm_i].rate;
        // Rate from reg.start to `start_pos - 1`:
        part_rate = partial_gamma_rate___(start_pos - 1, regions[gamm_i]);
        // Set `start_rate`:
        start_rate = cum_wt + part_rate;
    }

    /*
     Rate for ending position:
     */
    safe_get_mut(end_pos, mut_i);
    cum_wt += regions[gamm_i].rate;
    // Find the GammaRegion for the new end:
    while ((gamm_i < (regions.size() - 1)) && (regions[gamm_i].end < end_pos)) {
        gamm_i++;
        if (regions[gamm_i].rate > 0) cum_wt += regions[gamm_i].rate;
    }
    if (regions[gamm_i].rate < 0) stop("update should not end on deleted region (start)");
    cum_wt -= regions[gamm_i].rate;
    // Rate from reg.start to `end_pos`:
    part_rate = partial_gamma_rate___(end_pos, regions[gamm_i]);
    // Set `start_rate`:
    end_rate = cum_wt + part_rate;

    start_end_set = true;

    return;
}



uint32 LocationSampler::sample(pcg64& eng,
                               const uint32& start,
                               const uint32& end) {

    update_start_end(start, end);

    long double u = runif_ab(eng, start_rate, end_rate);

    // Find the GammaRegion:
    uint32 i = 0;
    long double cum_wt = 0;
    for (; i < regions.size(); i++) {
        if (regions[i].rate > 0) cum_wt += regions[i].rate;
        if (cum_wt > u) break;
    }

    /*
     Find the location within the Gamma region:
     */
    uint32 pos;
    gamma_sample(pos, u, cum_wt, i);

    return pos;
}
uint32 LocationSampler::sample(pcg64& eng) const {

    long double u = runif_01(eng) * total_rate;

    // Find the GammaRegion:
    uint32 i = 0;
    long double cum_wt = 0;
    for (; i < regions.size(); i++) {
        if (regions[i].rate > 0) cum_wt += regions[i].rate;
        if (cum_wt > u) break;
    }

    /*
     Find the location within the Gamma region:
     */
    uint32 pos;
    gamma_sample(pos, u, cum_wt, i);

    return pos;
}



/*
 Sample within a gamma region:
 */
inline void LocationSampler::gamma_sample(uint32& pos,
                                          long double& u,
                                          long double& cum_wt,
                                          const uint32& gam_i) const {

    if (regions[gam_i].rate < 0) stop("Cannot sample from a deleted region");

    const GammaRegion& reg(regions[gam_i]);
    const uint32& start(reg.start);
    const uint32& end(reg.end);


    // Update numbers so that `u` points to a place inside this region:
    cum_wt -= reg.rate;
    u -= cum_wt;
    u /= reg.gamma;
    cum_wt = 0;
    pos = start;

    /*
     If there are no mutations or if `end` is before the first mutation,
     then we don't need to use the `mutations` field at all.
    */
    if (var_seq->mutations.empty() || (end < var_seq->mutations.front().new_pos)) {

        for (; pos <= end; pos++) {
            cum_wt += nt_rates[var_seq->ref_seq->nucleos[pos]];
            if (cum_wt > u) break;
        }
        return;

    }

    // Index to the first Mutation object not past `start` position:
    uint32 mut_i = var_seq->get_mut_(start);
    // Current position
    pos = start;

    /*
     If `start` is before the first mutation (resulting in
     `mut_i == var_seq->mutations.size()`),
     we must pick up any nucleotides before the first mutation.
     */
    if (mut_i == var_seq->mutations.size()) {
        mut_i = 0;
        for (; pos < var_seq->mutations[mut_i].new_pos; pos++) {
            cum_wt += nt_rates[var_seq->ref_seq->nucleos[pos]];
            if (cum_wt > u) return;
        }
    }

    /*
     Now, for each subsequent mutation except the last, add all nucleotides
     at or after its position but before the next one.
     I'm adding `pos <= end` inside all while-statement checks to make sure
     it doesn't keep going after we've reached `end`.
     */
    uint32 next_mut_i = mut_i + 1;
    while (pos <= end && next_mut_i < var_seq->mutations.size()) {
        while (pos <= end && pos < var_seq->mutations[next_mut_i].new_pos) {
            char c = var_seq->get_char_(pos, mut_i);
            cum_wt += nt_rates[c];
            if (cum_wt > u) return;
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos <= end && pos < var_seq->seq_size) {
        char c = var_seq->get_char_(pos, mut_i);
        cum_wt += nt_rates[c];
        if (cum_wt > u) return;
        ++pos;
    }

    return;

}



