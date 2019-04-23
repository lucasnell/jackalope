
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
`ind` is the index to the current region in the vector of regions.
`erase_inds` stores indices for region(s) to be erased if the deletion
entirely spans one or more region(s).
Adding to this variable will result in the current region being erased.
*/
void GammaRegion::deletion_adjust(const uint32& ind,
                                  std::vector<uint32>& erase_inds,
                                  const uint32& del_start,
                                  const uint32& del_end,
                                  const sint32& del_size) {

    // Total overlap
    if ((del_start <= start) && (del_end >= end)) {
        erase_inds.push_back(ind);
        return;
    }
    // Deletion is totally inside this region but doesn't entirely overlap it
    if ((del_start > start) && (del_end < end)) {
        end += del_size;
        return;
    }
    // No overlap and deletion starts after it
    if (del_start > end) return;
    // No overlap and deletion starts before it
    if (del_end < start) {
        start += del_size;
        end += del_size;
        return;
    }
    // Partial overlap at the start
    if ((del_end >= start) && (del_start < start)) {
        start = del_start;
        end += del_size;
        return;
    }
    // Partial overlap at the end
    if ((del_start <= end) && (del_end > end)) {
        end = del_start - 1;
    }
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

    seq_size += static_cast<double>(size_change);
    uint32 idx = get_gamma_idx(pos);


    /*
    Insertions
    */
    if (size_change > 0) {
        regions[idx].end += size_change;
        idx++;
        // update all following ranges:
        while (idx < regions.size()) {
            regions[idx].end += size_change;
            regions[idx].start += size_change;
            idx++;
        }
        return;
    }

    /*
    Deletions
    */
    const uint32& del_start(pos);
    uint32 del_end = pos;
    del_end -= (size_change + 1);

    // Iterate through and adjust all regions including and following the deletion:
    std::vector<uint32> erase_inds;
    while (idx < regions.size()) {
        regions[idx].deletion_adjust(idx, erase_inds, del_start, del_end,
                                     size_change);
        idx++;
    }

    /*
    If any regions need erasing, their indices will be stored in erase_inds.
    They will be consecutive indices, so we only need to access the front and back.
    */
    if (erase_inds.size() == 1) {
        regions.erase(regions.begin() + erase_inds.front());
    } else if (erase_inds.size() > 1) {
        regions.erase(regions.begin() + erase_inds.front(),
                      regions.begin() + erase_inds.back() + 1);
    }

    return;
}






// To return the overall rate for an entire sequence or part of it:

double LocationSampler::calc_rate(uint32 start,
                                  uint32 end,
                                  const bool& ranged) const {

    double out = 0;

    if (var_seq->size() == 0) return out;

    if (!ranged) {
        start = 0;
        end = var_seq->size() - 1;
    }

    if ((var_seq->size() - 1) != regions.back().end) {
        stop("gammas and var_seq sizes don't match inside LocationSampler");
    }

    /*
    If there are no mutations or if `end` is before the first mutation,
    then we don't need to use the `mutations` field at all.
    (I'm using separate statements to avoid calling `front()` on an empty deque.)
    */
    bool use_mutations = true;
    if (var_seq->mutations.empty()) {
        use_mutations = false;
        if ((var_seq->ref_seq->nucleos.size() - 1) != regions.back().end) {
            stop("gammas and var_seq ref sizes don't match inside LocationSampler");
        }
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




double LocationSampler::deletion_rate_change(const sint32& size_mod,
                                             const uint32& start) {

    uint32 n_del = std::abs(size_mod); // # bp deleted

    uint32 end = start + n_del - 1;
    if (end >= var_seq->size()) end = var_seq->size() - 1;

    std::string seq;
    seq.reserve(n_del);
    uint32 mut_ = var_seq->get_mut_(start);
    var_seq->set_seq_chunk(seq, start, end - start + 1, mut_);

    double r, out = 0;
    uint32 seq_i = 0;
    uint32 idx = get_gamma_idx(start);

    while (seq_i < seq.size()) {
        GammaRegion& reg(regions[idx]);
        while (((seq_i + start) <= reg.end) && (seq_i < seq.size())) {
            r = reg.gamma * nt_rates[seq[seq_i]];
            out -= r;
            reg.rate -= r;
            total_rate -= r;
            seq_i++;
        }
        idx++;
    }

    return out;
}




uint32 LocationSampler::sample(pcg64& eng, const uint32& start, const uint32& end) {
    uint32 pos = (runif_01(eng) * (end - start + 1)) + start;
    return pos;
}
uint32 LocationSampler::sample(pcg64& eng) {
    uint32 pos = runif_01(eng) * var_seq->size();
    return pos;
}
