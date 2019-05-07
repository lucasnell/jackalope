
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
 Adjust positions (NOT rate, except for a full deletion) for a deletion.
 If the deletion totally overlaps it, this function changes this regions's rate
 to -1, and makes both start and end 0;
 it lastly returns true.
 Otherwise, it returns false.
*/
bool GammaRegion::deletion_adjust(const uint32& del_start,
                                  const uint32& del_end,
                                  const uint32& del_size,
                                  const uint32& i,
                                  RegionRejSampler& region_sampler) {

    if (deleted_) return false;

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
        // rate = -1;
        Rcout << "<del true> ";
        deleted_ = true;
        std::deque<uint32>& inds(region_sampler.levels[this->level].inds);
        auto iter = std::find(inds.begin(), inds.end(), i);
        if (iter == inds.end()) stop("i cannot be removed bc it's not here");
        inds.erase(iter);
        // for (uint32 j = 0; j < inds.size(); j++) if (inds[j] == i) stop("inds still has i");
        if (inds.size() == 0) {
            region_sampler.all_lvl_rate -= region_sampler.levels[this->level].lvl_rate;
            region_sampler.levels[this->level].lvl_rate = 0;
        }
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


void GammaRegion::check_rate(const double& threshold) {
    if (rate <= threshold) {
        Rcout << "<sub true> ";
        // rate = -1000;
        deleted_ = true;
    }
    return;
}







/*
 =======================================================================================
 =======================================================================================

 RejLevel and RegionRejSampler methods

 =======================================================================================
 =======================================================================================
 */




uint32 RejLevel::sample(pcg64& eng, const std::vector<GammaRegion>& regions) const {

    if (inds.size() == 0) stop("inds.size() == 0");

    double u;
    uint32 j;
    uint32 iters = 0;
    while (iters < 10) {
        j = runif_01(eng) * inds.size();
        u = runif_01(eng) * max_rate;
        if (u <= regions[inds[j]].rate) break;
        iters++;
    }
    return inds[j];
}
/*
 Add a Gamma region to this level.Used in constructing the object.
 A check is used before this function that keeps it from being used if
 the region's rate is <= 0 (thus is an invariant or deleted region), so no check for
 that is needed here.
 */
void RejLevel::add(std::vector<GammaRegion>& regions,
                   const uint32& i,
                   const uint32& lvl) {
    if (regions[i].rate > max_rate) max_rate = regions[i].rate;
    inds.push_back(i);
    lvl_rate += regions[i].rate;
    regions[i].level = lvl;
    return;
}
/*
 Remove a Gamma region from this level.
 Used in `LocationSampler::update_gamma_regions` when a deletion totally
 removes a Gamma region.
 It's not necessary for this function to update `*_rate` fields bc that's done
 inside the `deletion_rate_change` method.
 */
void RejLevel::rm(const uint32& i) {
    auto iter = std::find(inds.begin(), inds.end(), i);
    if (iter == inds.end()) stop("i cannot be removed bc it's not here");
    inds.erase(iter);
    return;
}

inline void RejLevel::reset_max(const std::vector<GammaRegion>& regions) {
    if (inds.size() == 0) {
        max_rate = 0;
        return;
    }
    max_rate = regions[inds[0]].rate;
    for (uint32 i = 1; i < inds.size(); i++) {
        if (regions[inds[i]].rate > max_rate) max_rate = regions[inds[i]].rate;
    }
    return;
}




RegionRejSampler::RegionRejSampler(std::vector<GammaRegion>& regions)
    : levels(0), all_lvl_rate(0) {

    if (regions.size() == 0) stop("cannot have 0-sized regions");
    if (regions.size() == 1) {

        all_lvl_rate = regions[0].rate;
        levels = std::vector<RejLevel>(1, RejLevel());
        levels[0].add(regions, 0, 0);

    } else {


        /*
         Go through once to get min and max rates:
         */
        double max_rate = 0;
        // (below, chose very large number that I'm sure is larger
        // than any/most rates)
        double min_rate = 1e10;
        for (uint32 i = 0; i < regions.size(); i++) {
            const GammaRegion& reg(regions[i]);
            if (reg.deleted()) continue; // <-- refers to deleted/invariant region
            if (reg.rate > max_rate) max_rate = reg.rate;
            if (reg.rate < min_rate) min_rate = reg.rate;
        }
        // Just in case...
        if (min_rate == 1e10) {
            stop("Min. rate that's > 0 is also > 1e10. This is entirely too high.");
        }

        /*
         Calculate how many levels we'll need, then add to `levels`:
         */
        uint32 n_levels = static_cast<uint32>(std::log2l(max_rate / min_rate)) + 1U;
        levels.reserve(n_levels);
        for (uint32 i = 0; i < n_levels; i++) levels.push_back(RejLevel());

        /*
         Now go back through regions and add each's info to proper level:
         */
        uint32 lvl_i;
        for (uint32 i = 0; i < regions.size(); i++) {
            // below refers to deleted/invariant region:
            if (regions[i].deleted()) continue;
            // If not one of those, then add it:
            lvl_i = static_cast<uint32>(std::log2l(regions[i].rate / min_rate));
            // lvl_i = 0;
            levels[lvl_i].add(regions, i, lvl_i);
            all_lvl_rate += regions[i].rate;
        }


    }

}




// Sample a Gamma region:
uint32 RegionRejSampler::sample(pcg64& eng,
                                const std::vector<GammaRegion>& regions) const {
    // Sample for the sampling-weight level:
    double u = runif_01(eng) * all_lvl_rate;
    double r = 0;
    uint32 i = 0;
    while (i < levels.size()) {
        r += levels[i].lvl_rate;
        if (u < r) break;
        i++;
    }
    if (i >= levels.size()) {
        Rcout << std::endl << std::endl << i << ' ' << levels.size() << std::endl;
        Rcout << r << ' ' << u << ' ' << all_lvl_rate << std::endl << std::endl;
        stop("i >= levels.size() in rej. region sampling");
    }
    // Sample for Gamma region within that level:
    uint32 pos = levels[i].sample(eng, regions);
    return pos;
};







/*
 ========================================================================================
 ========================================================================================

 LocationSampler methods

 ========================================================================================
 ========================================================================================
 */




/*
 Based on a sequence position, return an index to the Gamma region it's inside.
 */
inline uint32 LocationSampler::get_gamma_idx(const uint32& pos) const {
    uint32 idx = pos * (static_cast<double>(regions.size()) /
        static_cast<double>(var_seq->size()));
    if (idx >= regions.size()) idx = regions.size() - 1;
    while (regions[idx].end < pos || regions[idx].deleted()) idx++;
    while (regions[idx].start > pos || regions[idx].deleted()) idx--;
    return idx;
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

    double gamma = gamma_mat(i,1); // this will be same even if it gets split

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

    /*
     Vector of # bases per GammaRegion, if splitting is needed for non-rejection
     sampling:
     */
    if (!rej_sample && n_regs > 1) {
        sizes = split_int(size, n_regs);
    } else {
        sizes = std::vector<uint32>(1, size);
    }

    // String to store sequence info:
    std::string seq;
    seq.reserve(sizes[0] + 1);

    double rate;

    for (uint32 j = 0; j < sizes.size(); j++) {

        // Calculate rate:
        rate = 0;
        // (in `set_seq_chunk` below, `seq` gets cleared before it's filled)
        var_seq->set_seq_chunk(seq, start, sizes[j], mut_i);
        for (const char& c : seq) rate += nt_rates[c];
        rate *= gamma;

        // Set LocationSampler info:
        total_rate += rate;
        regions.push_back(GammaRegion(gamma, start, start + sizes[j] - 1, rate));

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
    if (rej_sample) {
        regions.reserve(gamma_mat.n_rows);
    } else regions.reserve((var_seq->size() / gamma_size) * 1.5);


    uint32 mut_i = 0;
    std::vector<uint32> sizes;
    if (!rej_sample) sizes.reserve(1000); // arbitrarily chosen

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
            if (!regions[idx].deleted()) {
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
    // uint32 idx = pos * (static_cast<double>(regions.size()) /
    //     static_cast<double>(var_seq->size()));
    // if (idx >= regions.size()) idx = regions.size() - 1;
    // while (regions[idx].end < pos || regions[idx].deleted()) idx++;
    // while (regions[idx].start > pos || regions[idx].deleted()) idx--;

    const uint32& del_start(pos);
    uint32 del_size = std::abs(size_change);
    uint32 del_end = pos + del_size - 1;

    // Iterate through and adjust all regions including and following the deletion:
    // bool weird = false;
    bool rm;
    for (; idx < regions.size(); idx++) {
        if (regions[idx].end < pos || regions[idx].deleted()) continue;
        // Next line returns true if the region is totally deleted:
        rm = regions[idx].deletion_adjust(del_start, del_end, del_size,
                                          idx, region_sampler);
        if (rm) {
            std::deque<uint32>& inds(region_sampler.levels[regions[idx].level].inds);
            for (uint32 j = 0; j < inds.size(); j++) {
                if (inds[j] == idx) stop("inds still has idx");
            }
        }
        // weird = regions[idx].rate <= 0 && !regions[idx].deleted();
        // if (rm || weird) {
        if (!rm && regions[idx].deleted())  stop("wtf");
        // if (rm || regions[idx].deleted()) {
        //     // If it's deleted, then remove it from sampling:
        //     RejLevel& lvl(region_sampler.levels[regions[idx].level]);
        //     lvl.rm(idx);
        // }
    }

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

double LocationSampler::calc_rate__(const uint32& start, const uint32& end) const {

    double out = 0;

    if (var_seq->size() == 0) return out;

    /*
     If there are no mutations or if `end` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     */
    if (var_seq->mutations.empty() || (end < var_seq->mutations.front().new_pos)) {

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
    while (pos <= end && pos < var_seq->seq_size) {
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






double LocationSampler::substitution_rate_change(const char& c, const uint32& pos) {
    uint32 i = get_gamma_idx(pos);
    GammaRegion& reg(regions[i]);
    RejLevel& lvl(region_sampler.levels[reg.level]);
    /*
     If the region affected by the substitution was the max at the sampling level
     and its rate is reduced by the substitution,
     then we must re-check that level to find the new max.
     This shouldn't happen too often.
     */
    bool was_max = (lvl.max_rate <= reg.rate);

    char c0 = var_seq->get_nt(pos);
    double d_rate = nt_rates[c] - nt_rates[c0];

    if (d_rate != 0) {

        d_rate *= reg.gamma;

        // Add to the rate for `reg` and check to see if that makes it zero
        reg.rate += d_rate;
        reg.check_rate();

        total_rate += d_rate;
        end_rate += d_rate;
        // If it was the max but decreased in rate, check that it's still the max:
        if (was_max && d_rate < 0) lvl.reset_max(regions);
        /*
         Also check to make sure it isn't the new max if it increased in rate:
         */
        if (lvl.max_rate < reg.rate) lvl.max_rate = reg.rate;
        // Update total rate for that level and for all levels:
        lvl.lvl_rate += d_rate;
        region_sampler.all_lvl_rate += d_rate;
        // Now remove from sampler if this mutation caused it to have a rate of 0
        if (reg.deleted()) {
            Rcout << "< sub rm > ";
            lvl.rm(i);
            if (lvl.inds.size() == 0 && lvl.lvl_rate != 0) {
                Rcout << "weird " << lvl.lvl_rate << ' ';
                region_sampler.all_lvl_rate -= lvl.lvl_rate;
                lvl.lvl_rate = 0;
            }
        }

    }

    return d_rate;
}

double LocationSampler::insertion_rate_change(const std::string& seq, const uint32& pos) {
    GammaRegion& reg(regions[get_gamma_idx(pos)]);
    RejLevel& lvl(region_sampler.levels[reg.level]);
    double gamma = reg.gamma;
    double d_rate = 0;
    for (const char& c : seq) d_rate += nt_rates[c];
    if (d_rate != 0) {
        d_rate *= gamma;
        reg.rate += d_rate;
        total_rate += d_rate;
        end_rate += d_rate;
        /*
         Because insertions will only ever increase the region's rate, we only need
         to check to make sure it isn't the new max:
         */
        if (lvl.max_rate < reg.rate) lvl.max_rate = reg.rate;
        // Update total rate for that level and for all levels:
        lvl.lvl_rate += d_rate;
        region_sampler.all_lvl_rate += d_rate;
    }
    return d_rate;
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

    double r, d_rate = 0;
    uint32 seq_i = 0;
    uint32 idx = get_gamma_idx(start);

    while (seq_i < seq.size()) {

        // Skip over deleted/invariant region if necessary:
        if (regions[idx].deleted()) {
            idx++;
            continue;
        }

        GammaRegion& reg(regions[idx]);
        RejLevel& lvl(region_sampler.levels[reg.level]);
        /*
         If any region affected by the deletion was the max at the sampling level,
         then we must re-check that level to find the new max.
         This shouldn't happen too often.
         */
        bool recheck_lvl_max = (lvl.max_rate == reg.rate);
        r = 0;
        while (((seq_i + start) <= reg.end) && (seq_i < seq.size())) {
            r += nt_rates[seq[seq_i]];
            seq_i++;
        }
        r *= (-1.0 * reg.gamma);
        d_rate += r;

        // Add to the rate for `reg` and check to see if that makes it zero
        reg.rate += d_rate;
        reg.check_rate();

        total_rate += r;
        end_rate += r;
        lvl.lvl_rate += r;
        region_sampler.all_lvl_rate += r;
        if (recheck_lvl_max) lvl.reset_max(regions);
        // If this deletion totally removed this region, remove it from sampling:
        if (reg.deleted()) {
            Rcout << "< del rm > ";
            lvl.rm(idx);
            for (uint32 ii = 0; ii < region_sampler.levels.size(); ii++) {
                RejLevel& lvl_(region_sampler.levels[ii]);
                auto iii = std::find(lvl_.inds.begin(), lvl_.inds.end(), idx);
                if (iii != lvl_.inds.end()) stop("duplicate idx found");
            }
            if (lvl.inds.size() == 0 && lvl.lvl_rate != 0) {
                Rcout << std::endl << std::endl << region_sampler.all_lvl_rate <<
                    ' ' << lvl.lvl_rate << ' ';
                double oo = 0;
                for (const RejLevel& l : region_sampler.levels) oo += l.lvl_rate;
                Rcout << oo << std::endl;
                stop("del rm making sized-0 lvl inds not have rate of 0");
                region_sampler.all_lvl_rate -= lvl.lvl_rate;
                lvl.lvl_rate = 0;
            }
        }
        idx++;
    }


    return d_rate;
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
inline double LocationSampler::partial_gamma_rate___(
        const uint32& end,
        const GammaRegion& reg) const {

    double out = 0;

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

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    stop("start end doesn't work");

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

    double cum_wt = 0;
    uint32 gamm_i = 0;
    // Skip over deleted region(s) if necessary:
    while (regions[gamm_i].deleted()) gamm_i++;
    double part_rate = 0;

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

    double u = runif_ab(eng, start_rate, end_rate);

    // Find the GammaRegion:
    uint32 i = 0;
    double cum_wt = 0;
    for (; i < regions.size(); i++) {
        if (regions[i].rate > 0) cum_wt += regions[i].rate;
        if (cum_wt > u) break;
    }

    /*
     Find the location within the Gamma region:
     */
    uint32 pos;
    cdf_region_sample(pos, u, cum_wt, i);
    // if (rej_sample) {
    //     rej_region_sample(pos, start, end, eng, i);
    // } else {
    //     cdf_region_sample(pos, u, cum_wt, i);
    // }


    // uint32 pos;
    //
    // if (rej_sample) {
    //
    //     // uint32 i = runif_01(eng) * regions.size();
    //     // while (regions[i].rate <= 0) i = runif_01(eng) * regions.size();
    //
    //     uint32 i = region_sampler.sample(eng);
    //
    //     rej_region_sample(pos, eng, i);
    //
    //     if (pos >= var_seq->size()) {
    //         Rcout << std::endl << pos << ' ' << var_seq->size() << std::endl;
    //         stop("pos > var_seq->size()");
    //     }
    //
    // } else {
    //
    //     double u = runif_01(eng) * total_rate;
    //
    //     // Find the GammaRegion:
    //     uint32 i = 0;
    //     double cum_wt = 0;
    //     for (; i < regions.size(); i++) {
    //         if (regions[i].rate > 0) cum_wt += regions[i].rate;
    //         if (cum_wt > u) break;
    //     }
    //
    //     cdf_region_sample(pos, u, cum_wt, i);
    //
    //     if (pos >= var_seq->size()) {
    //         Rcout << std::endl << pos << ' ' << var_seq->size() << std::endl;
    //         stop("pos > var_seq->size() in cdf");
    //     }
    //
    // }

    return pos;
}
uint32 LocationSampler::sample(pcg64& eng) const {

    uint32 pos;

    // for (uint32 i = 0; i < region_sampler.levels.size(); i++) {
    //     const std::deque<uint32>& inds(region_sampler.levels[i].inds);
    //     for (uint32 j = 0; j < inds.size(); j++) {
    //         uint32 k = inds[j];
    //         if (regions[k].deleted() || regions[k].rate <= 0) {
    //             Rcout << std::endl << regions[k].rate << ' ' << regions[k].deleted() << std::endl;
    //             stop("regions[k].rate <= 0 || deleted in LocationSampler::sample");
    //         }
    //     }
    // }


    uint32 i = region_sampler.sample(eng, regions);

    // if (regions[i].rate <= 0 || regions[i].deleted()) {
    //     Rcout << std::endl << regions[i].rate << ' ' << regions[i].deleted() << std::endl;
    //     stop("regions[i] <= 0");
    // }
    // if (i > regions.size()) stop("i > regions.size()");

    if (rej_sample) {

        rej_region_sample(pos, eng, i);

        // if (pos >= var_seq->size()) {
        //     Rcout << std::endl << pos << ' ' << var_seq->size() << std::endl;
        //     stop("pos > var_seq->size()");
        // }

    } else {

        cdf_region_sample(pos, eng, i);

        // if (pos >= var_seq->size()) {
        //     Rcout << std::endl << pos << ' ' << var_seq->size() << std::endl;
        //     stop("pos > var_seq->size() in cdf");
        // }

    }

    return pos;
}



/*
 Sample within one GammaRegion using rejection method:
 */
inline void LocationSampler::rej_region_sample(uint32& pos,
                                               pcg64& eng,
                                               const uint32& gam_i) const {

    if (regions[gam_i].rate <= 0) stop("Cannot sample from a deleted/invariant region");

    const GammaRegion& reg(regions[gam_i]);
    // Make refs for `start` and `end` for this region
    const uint32& start(reg.start);
    const uint32& end(reg.end);

    // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // if (start > end) {
    //     Rcout << std::endl << start << ' ' << end << std::endl;
    //     stop("start > end");
    // }
    // // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (start == end) {
        pos = start;
        return;
    }

    double u, q;
    char c;
    uint32 size_ = end + 1;
    size_ -= start;
    uint32 n_iters = 0;

    while (n_iters < 100) {
        pos = runif_01(eng) * size_;
        pos += start;
        // // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // if (pos > var_seq->size()) {
        //     Rcout << std::endl << pos << ' ' << start <<
        //         ' ' << end << ' ' << size_ << ' ' << var_seq->size() << std::endl;
        //     stop("pos > var_seq->size() #1");
        // }
        // // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        c = var_seq->get_nt(pos);
        q = nt_rates[c];
        u = runif_01(eng) * max_q;
        if (u <= q) break;
        n_iters++;
    }

    return;

}
/*
 Sample within one GammaRegion using cdf (non-rejection) method:
 */
inline void LocationSampler::cdf_region_sample(uint32& pos,
                                               double& u,
                                               double& cum_wt,
                                               const uint32& gam_i) const {

    if (regions[gam_i].rate <= 0) stop("Cannot sample from a deleted region");

    const GammaRegion& reg(regions[gam_i]);
    const uint32& start(reg.start);
    const uint32& end(reg.end);

    if (start == end) {
        pos = start;
        return;
    }

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


/* Same as above, but when using rejection sampling for region */
inline void LocationSampler::cdf_region_sample(uint32& pos,
                                               pcg64& eng,
                                               const uint32& gam_i) const {

    if (regions[gam_i].rate <= 0) stop("Cannot sample from a deleted region");

    const GammaRegion& reg(regions[gam_i]);
    const uint32& start(reg.start);
    const uint32& end(reg.end);

    if (end >= var_seq->size() || start >= var_seq->size()) {
        Rcout << std::endl << start << ' ' << end << ' ' << var_seq->size() << std::endl;
        stop("end or start >= var_seq->size() inside cdf");
    }

    if (start == end) {
        pos = start;
        return;
    }

    // Sample for a place inside this region:
    double cum_wt = 0;
    double u = runif_01(eng) * (reg.rate / reg.gamma);
    pos = start;

    /*
     If there are no mutations or if `end` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     */
    if (var_seq->mutations.empty() || (end < var_seq->mutations.front().new_pos)) {

        for (; pos <= end; pos++) {
            cum_wt += nt_rates[var_seq->ref_seq->nucleos[pos]];
            if (cum_wt > u) {
                if (pos > end) {
                    Rcout << std::endl << pos << ' ' << end << std::endl;
                    stop("pos > end in cdf - #0");
                }
                break;
            }
        }
        if (pos > end) {
            Rcout << std::endl << pos << ' ' << end << std::endl;
            stop("pos > end in cdf - #1");
        }
        return;

    }

    // Index to the first Mutation object not past `start` position:
    uint32 mut_i = var_seq->get_mut_(start);

    /*
     If `start` is before the first mutation (resulting in
     `mut_i == var_seq->mutations.size()`),
     we must pick up any nucleotides before the first mutation.
     */
    if (mut_i == var_seq->mutations.size()) {
        mut_i = 0;
        for (; pos < var_seq->mutations[mut_i].new_pos; pos++) {
            cum_wt += nt_rates[var_seq->ref_seq->nucleos[pos]];
            if (cum_wt > u) {
                if (pos > end) {
                    Rcout << std::endl << pos << ' ' << end << std::endl;
                    stop("pos > end in cdf - #2");
                }
                return;
            }
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
            if (cum_wt > u) {
                if (pos > end) {
                    Rcout << std::endl << pos << ' ' << end << std::endl;
                    stop("pos > end in cdf - #3");
                }
                return;
            }
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos <= end && pos < var_seq->seq_size) {
        char c = var_seq->get_char_(pos, mut_i);
        cum_wt += nt_rates[c];
        if (cum_wt > u) {
            if (pos > end) {
                Rcout << std::endl << pos << ' ' << end << std::endl;
                stop("pos > end in cdf - #4");
            }
            return;
        }
        ++pos;
    }

    // if (pos > end) {
    //     Rcout << std::endl << pos << ' ' << end << std::endl;
    //     Rcout << cum_wt << ' ' << (reg.rate / reg.gamma) << std::endl;
    //     stop("pos > end in cdf - end");
    // }

    return;

}


