
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
 Adjust bounds for a deletion.
 It changes the `deleted` field to true if this region should be deleted.
 */
void Region::del_adjust_bounds(const uint64& del_start,
                               const uint64& del_end) {

    // No overlap and deletion starts after it
    if (del_start > end) return;

    uint64 del_size = del_end - del_start + 1;

    // No overlap and deletion starts before it
    if (del_end < start) {
        start -= del_size;
        end -= del_size;
        return;
    }
    /*
     Total overlap
     */
    if ((del_start <= start) && (del_end >= end)) {
        start = del_start - 1;
        end = del_start - 1;
        deleted = true;
        return;
    }
    /*
     Deletion is totally inside this region but doesn't entirely overlap it
     or touch any of the end points
     */
    if ((del_start > start) && (del_end < end)) {
        end -= del_size;
        return;
    }
    // Partial overlap at the start
    if ((del_end >= start) && (del_start <= start)) {
        start = del_start;
        end -= del_size;
        return;
    }
    // Partial overlap at the end
    if ((del_start <= end) && (del_end >= end)) {
        end = del_start - 1;
    }
    return;
}




inline void RegionTree::construct_tips_one_row(const arma::mat& gamma_mat,
                                               const uint64& region_size,
                                               const VarSequence* var_seq,
                                               const std::vector<double>& nt_rates,
                                               const uint64& i,
                                               uint64& mut_i,
                                               std::vector<uint64>& sizes) {

    double gamma = gamma_mat(i,1); // this will be same even if it gets split
    if (gamma <= 0) return; // no need to add invariant regions

    uint64 start, size, n_regs;
    if (i > 0) {
        /*
         Even though `gamma_mat.col(0)` uses 1-based indices, I'm not subtracting 1
         below because I'm looking at previous row:
         */
        start = static_cast<uint64>(gamma_mat(i-1,0));
    } else start = 0;

    /*
     Because `size = end - start + 1` and bc `end = gamma_mat(i,0) - 1`
     (`gamma_mat.col(0)` uses 1-based indices),
     we don't need to add or subtract 1.
     */
    size = static_cast<uint64>(gamma_mat(i,0)) - start;

    n_regs = std::round(static_cast<double>(size) / static_cast<double>(region_size));

    // Vector of # bases per Region
    if (n_regs > 1) {
        sizes = split_int(size, n_regs);
    } else {
        sizes = std::vector<uint64>(1, size);
    }

    // String to store sequence info:
    std::string seq;
    seq.reserve(sizes[0] + 1);
    double rate;

    for (uint64 j = 0; j < sizes.size(); j++) {

        // Calculate rate:
        rate = 0;
        // (in `set_seq_chunk` below, `seq` gets cleared before it's filled)
        var_seq->set_seq_chunk(seq, start, sizes[j], mut_i);
        for (const char& c : seq) rate += nt_rates[c];
        rate *= gamma;

        // Set RegionTree info:
        total_rate += rate;
        tips.push_back(Region(gamma, start, start + sizes[j] - 1, rate));

        // Update start for next iteration:
        start += sizes[j];
    }


    return;

}





void RegionTree::construct_tips(arma::mat gamma_mat,
                                const uint64& region_size,
                                const VarSequence* var_seq,
                                const std::vector<double>& nt_rates) {

    total_rate = 0;

    tips.clear();
    tips.reserve((var_seq->size() / region_size) * 1.5);

    // Sort from first to last region
    arma::uvec sort_inds = arma::sort_index(gamma_mat.col(0));
    gamma_mat = gamma_mat.rows(sort_inds);

    uint64 mut_i = 0;
    std::vector<uint64> sizes;
    sizes.reserve(1000); // arbitrarily chosen

    for (uint64 i = 0; i < gamma_mat.n_rows; i++) {
        construct_tips_one_row(gamma_mat, region_size, var_seq, nt_rates, i, mut_i, sizes);
    }

    clear_memory<std::vector<Region>>(tips);  // in case I reserved too much memory

    return;
}



void RegionTree::construct_nodes() {

    uint64 n_tips = tips.size();
    if (n_tips == 0) stop("no tips present");
    // If just one tip, then we can keep this at size 0:
    if (n_tips == 1) {
        nodes.clear();
        return;
    }

    // Else, let's figure out how many levels we need first:
    uint64 n_lvls = std::ceil(std::log2(static_cast<double>(n_tips)));

    // Resize `nodes` and all inside levels by iterating backwards through each level:
    nodes.resize(n_lvls);
    double n = n_tips;
    for (uint64 i = 1; i <= n_lvls; i++) {
        n = std::ceil(n / 2.0);
        // Using this type of construction to make sure they're all zeroes
        nodes[n_lvls-i] = std::vector<double>(static_cast<uint64>(n), 0.0);
    }
    // Fill all node values:
    for (uint64 tip_i = 0; tip_i < n_tips; tip_i++) {
        double rate = tips[tip_i].rate;
        // update nodes from a tip:
        update_nodes(tip_i, rate);
    }

    return;

}



const Region* RegionTree::search(double& u) const {

    if (nodes.size() == 0) return &tips[0];

    // Traverse down tree to see which tip to return:
    uint64 node_i = 0;
    bool go_right;
    // Go up tree and update all necessary nodes:
    for (uint64 lvl_i = 0; lvl_i < nodes.size(); lvl_i++) {
        const double& threshold(nodes[lvl_i][node_i]);
        go_right = u > threshold;
        node_i <<= 1ULL; // traversing down means using this type of bit-shift
        if (go_right) {
            u -= threshold;
            node_i++;
        }
    }

    mut_tip_ = node_i; // saving this to later update it for the mutation

    return &tips[node_i];

}





void LocationSampler::update(const double& d_rate,
                             const sint64& size_change,
                             const uint64& pos) {

    // Update `end_rate` and `end_pos` in `bounds`
    bounds.update(d_rate, size_change);

    // Substitutions
    if (size_change == 0) {
        regions.sub_update(d_rate);
        return;
    }

    // Insertions
    if (size_change > 0) {
        regions.ins_update(d_rate, static_cast<uint64>(size_change));
        return;
    }

    // Deletions
    uint64 end = pos - size_change - 1;
    regions.del_update(d_rate, pos, end, del_rate_changes);

    return;
}











double LocationSampler::deletion_rate_change(const uint64& del_size,
                                             const uint64& start) const {

    uint64 end = start + del_size - 1;
    if (end >= var_seq->size()) end = var_seq->size() - 1;

    // Prep `del_rate_changes` to store rate changes for deletions spanning >1 regions:
    del_rate_changes.clear();

    std::string seq;
    seq.reserve(del_size);
    uint64 mut_i;
    safe_get_mut(start, mut_i);

    var_seq->set_seq_chunk(seq, start, end - start + 1, mut_i);

    double r, out = 0;
    uint64 seq_i = 0;
    uint64 idx = regions.mut_tip();

    while (seq_i < seq.size()) {
        const Region& reg(regions.tips[idx]);
        if (!reg.deleted) {
            r = 0;
            while (((seq_i + start) <= reg.end) && (seq_i < seq.size())) {
                r -= reg.gamma * nt_rates[seq[seq_i]];
                seq_i++;
            }
            out += r;
            del_rate_changes.push_back(r);
        } else del_rate_changes.push_back(0);
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
inline void LocationSampler::safe_get_mut(const uint64& pos, uint64& mut_i) const {

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
 It just returns the rate (inclusive) from `reg.start` to `end`
 */
inline double LocationSampler::partial_gamma_rate___(
        const uint64& end,
        const Region& reg) const {

    double out = 0;

    if (end > reg.end) stop("end > reg.end");
    if (end < reg.start) stop("end < reg.start");
    if (end == reg.end) return reg.rate;

    /*
     If there are no mutations or if `end` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     */
    if (var_seq->mutations.empty() || (end < var_seq->mutations.front().new_pos)) {

        for (uint64 i = reg.start; i <= end; i++) {
            out += nt_rates[var_seq->ref_seq->nucleos[i]];
        }
        out *= reg.gamma;
        return out;

    }

    // Current position
    uint64 pos = reg.start;
    // Index to the first Mutation object not past `pos`:
    uint64 mut_i = var_seq->get_mut_(pos);

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
    uint64 next_mut_i = mut_i + 1;
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
void LocationSampler::new_bounds(const uint64& start,
                                 const uint64& end) {

    uint64& start_pos(bounds.start_pos);
    uint64& end_pos(bounds.end_pos);
    double& start_rate(bounds.start_rate);
    double& end_rate(bounds.end_rate);
    bool& start_end_set(bounds.start_end_set);

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
        end_rate = regions.total_rate;
        return;
    }

    double cum_wt = 0;
    uint64 gamm_i = 0;
    double part_rate = 0;

    uint64 mut_i = 0;

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
        gamm_i = 0;
        cum_wt = regions.tips[gamm_i].rate;
        // Find the Region for the new start:
        while ((gamm_i < (regions.tips.size() - 1)) &&
               (regions.tips[gamm_i].end < (start_pos - 1))) {
            gamm_i++;
            cum_wt += regions.tips[gamm_i].rate;
        }
        cum_wt -= regions.tips[gamm_i].rate;
        // Rate from reg.start to `start_pos - 1`:
        part_rate = partial_gamma_rate___(start_pos - 1, regions.tips[gamm_i]);
        // Set `start_rate`:
        start_rate = cum_wt + part_rate;
    }

    /*
     Rate for ending position:
     */
    safe_get_mut(end_pos, mut_i);
    cum_wt += regions.tips[gamm_i].rate;
    // Find the Region for the new end:
    while ((gamm_i < (regions.tips.size() - 1)) && (regions.tips[gamm_i].end < end_pos)) {
        gamm_i++;
        cum_wt += regions.tips[gamm_i].rate;
    }
    cum_wt -= regions.tips[gamm_i].rate;
    // Rate from reg.start to `end_pos`:
    part_rate = partial_gamma_rate___(end_pos, regions.tips[gamm_i]);
    // Set `start_rate`:
    end_rate = cum_wt + part_rate;

    start_end_set = true;

    return;
}



uint64 LocationSampler::sample(pcg64& eng,
                               const uint64& start,
                               const uint64& end) {

    new_bounds(start, end);

    double u = runif_ab(eng, bounds.start_rate, bounds.end_rate);

    // Find the Region:
    const Region* reg = regions.search(u);

    /*
     Find the location within the Gamma region:
     */
    uint64 pos;
    cdf_region_sample(pos, u, reg);

    return pos;
}
uint64 LocationSampler::sample(pcg64& eng) const {

    double u = runif_01(eng) * regions.total_rate;

    // Find the Region:
    const Region* reg = regions.search(u);

    /*
     Find the location within the Gamma region:
     */
    uint64 pos;
    cdf_region_sample(pos, u, reg);

    return pos;
}



/*
 Sample within a gamma region using CDF method:
 */
inline void LocationSampler::cdf_region_sample(uint64& pos,
                                               double& u,
                                               const Region* reg) const {

    const uint64& start(reg->start);
    const uint64& end(reg->end);

    /*
     Because `reg->rate` incorporates `reg->gamma`, we'd have to multiply all our
     nucleotide-level rates by `reg->gamma` if we didn't run the following line:
     */
    u /= reg->gamma;

    pos = start;
    double cum_wt = 0;

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
    uint64 mut_i = var_seq->get_mut_(start);
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
    uint64 next_mut_i = mut_i + 1;
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



