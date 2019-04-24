
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
                                  const uint32& del_size) {

    // No overlap and deletion starts after it
    if (del_start > end) return;
    // No overlap and deletion starts before it
    if (del_end < start) {
        start -= del_size;
        end -= del_size;
        return;
    }
    // Total overlap
    if ((del_start <= start) && (del_end >= end)) {
        erase_inds.push_back(ind);
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





void LocationSampler::construct_gammas(arma::mat gamma_mat) {

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
    regions.reserve(gamma_mat.n_rows);

    double gamma;
    uint32 start;
    uint32 end;
    double rate;
    uint32 mut_i = 0;

    for (uint32 i = 0; i < gamma_mat.n_rows; i++) {

        // Below, I'm subtracting 1 to go back to 0-based indexing
        end = static_cast<uint32>(gamma_mat(i,0)) - 1;

        if (i > 0) {
            start = static_cast<uint32>(gamma_mat(i-1,0));  // not subtracting -1
        } else start = 0;

        gamma = gamma_mat(i,1);

        std::string seq;
        seq.reserve(end - start + 1);
        rate = 0;
        var_seq->set_seq_chunk(seq, start, end - start + 1, mut_i);
        for (const char& c : seq) rate += nt_rates[c];

        regions.push_back(GammaRegion(gamma, start, end, rate));

        total_rate += rate;
    }

    if (regions.back().end >= var_seq->size()) {
        stop("on setup, regions.back().end >= var_seq->size()");
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

        if (regions.back().end != (var_seq->size() - 1)) {
            Rcout << std::endl <<  "INS (region, seq) " << regions.back().end << ", " <<
                var_seq->size() << std::endl;
            stop("while running, regions.back().end != (var_seq->size() - 1)");
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
    std::vector<uint32> erase_inds;
    while (idx < regions.size()) {
        regions[idx].deletion_adjust(idx, erase_inds, del_start, del_end,
                                     del_size);
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


    if (regions.back().end != (var_seq->size() - 1)) {
        Rcout << std::endl <<  "DEL (region, seq) " << regions.back().end << ", " <<
            var_seq->size() << std::endl;
        stop("while running, regions.back().end != (var_seq->size() - 1)");
    }

    return;
}






// To return the overall rate for an entire sequence or part of it:

// Inner method that does most of the work for `calc_rate` below

double LocationSampler::calc_rate__(uint32 start, uint32 end) const {

    double out = 0;

    if (var_seq->size() == 0) return out;

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








uint32 LocationSampler::sample(pcg64& eng,
                               const uint32& start,
                               const uint32& end) const {
    uint32 pos = (runif_01(eng) * (end - start + 1)) + start;
    return pos;
}
uint32 LocationSampler::sample(pcg64& eng) const {

    long double u = runif_01(eng) * total_rate;

    // Find the GammaRegion:
    uint32 i = 0;
    long double cum_wt = 0;
    for (; i < regions.size(); i++) {
        cum_wt += regions[i].rate;
        if (cum_wt > u) break;
    }

    /*
     Find the location within the Gamma region:
     */
    uint32 pos = gamma_sample(u, cum_wt, i);

    return pos;
}




inline uint32 LocationSampler::gamma_sample(long double& u,
                                            long double& cum_wt,
                                            const uint32& gam_i) const {

    const GammaRegion& reg(regions[gam_i]);

    // Fill sequence of Gamma region:
    std::string reg_seq;
    uint32 reg_size = reg.end - reg.start + 1;
    reg_seq.reserve(reg_size);
    uint32 mut_i;
    safe_get_mut(reg.start, mut_i);
    var_seq->set_seq_chunk(reg_seq, reg.start, reg_size, mut_i);

    // Now see where `u` points to:
    cum_wt -= reg.rate;
    u -= cum_wt;
    u /= reg.gamma;
    uint32 pos = 0;
    cum_wt = 0;
    for (; pos < reg_size; pos++) {
        cum_wt += nt_rates[reg_seq[pos]];
        if (cum_wt > u) break;
    }

    pos += reg.start;

    return pos;

}



