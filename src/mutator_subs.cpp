

#include "mutator_subs.h" // SubMutator and debugging preprocessor directives



#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class


#include "var_classes.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "util.h"  // str_stop





void SubMutator::new_chrom(VarChrom& var_chrom_, pcg64& eng) {

    var_chrom = &var_chrom_;

    // (Gammas go from 0 to (n-1), invariants are n.)
    const uint8 n = Q.size();

    const uint64 N = var_chrom_.size();
    const uint64 N0 = rate_inds.size();

    if (invariant <= 0) {

        if (N0 < N) {

            for (uint64 i = 0; i < N0; i++) {
                rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            }
            for (uint64 i = N0; i < N; i++) {
                rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
            }

        } else {

            if (N0 > N) rate_inds.resize(N);
            for (uint64 i = 0; i < N; i++) {
                rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            }

        }

    } else {

        if (N0 < N) {

            for (uint64 i = 0; i < N0; i++) {
                if (runif_01(eng) > invariant) {
                    rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
                } else rate_inds[i] = n;
            }
            for (uint64 i = N0; i < N; i++) {
                if (runif_01(eng) > invariant) {
                    rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
                } else rate_inds.push_back(n);
            }

        } else {

            if (N0 > N) rate_inds.resize(N);
            for (uint64 i = 0; i < N; i++) {
                if (runif_01(eng) > invariant) {
                    rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
                } else rate_inds[i] = n;
            }

        }
    }

    return;
}




void SubMutator::new_branch(const double& b_len) {

    // UNREST model
    if (U.size() == 0) {
        for (uint32 i = 0; i < Q.size(); i++) {
            // Adjust P(t) matrix
            Pt_calc(Q[i], 30, b_len, Pt[i]);
            // Now adjust the alias samplers:
            std::vector<AliasSampler>& samp(samplers[i]);
            for (uint32 j = 0; j < 4; j++) {
                samp[j] = AliasSampler(Pt[i].row(j));
            }
        }
    } else {
#ifdef __JACKALOPE_DEBUG
        if (U.size() != Q.size()) stop("SubMutator::new_branch-> U.size() != Q.size()");
        if (Ui.size() != Q.size()) stop("SubMutator::new_branch-> Ui.size() != Q.size()");
        if (L.size() != Q.size()) stop("SubMutator::new_branch-> L.size() != Q.size()");
#endif
        // All other models
        for (uint32 i = 0; i < Q.size(); i++) {
            // Adjust P(t) matrix
            Pt_calc(U[i], Ui[i], L[i], b_len, Pt[i]);
            // Now adjust the alias samplers:
            std::vector<AliasSampler>& samp(samplers[i]);
            for (uint32 j = 0; j < 4; j++) {
                samp[j] = AliasSampler(Pt[i].row(j));
            }
        }
    }

    return;

}

void SubMutator::add_subs(const double& b_len, pcg64& eng) {

    new_branch(b_len);

    uint32 max_gamma = Q.size() - 1;
    std::string bases = "TCAG";

    /*
     If there are no mutations, then we obviously don't need to use the `mutations` field.
     */
    if (var_chrom->mutations.empty()) {

        for (uint64 pos = 0; pos < var_chrom->size(); pos++) {

            uint8& rate_i(rate_inds[pos]);
            if (rate_i > max_gamma) continue; // this is an invariant region

            uint8 c_i = char_map[var_chrom->ref_chrom->nucleos[pos]];
            if (c_i > 3) continue; // only changing T, C, A, or G
            AliasSampler& samp(samplers[rate_i][c_i]);
            uint8 nt_i = samp.sample(eng);
            if (nt_i != c_i) var_chrom->add_substitution(bases[nt_i], pos);

        }

        return;

    }

    /*
     *********************************************************************************
     *********************************************************************************
     left off here
     *********************************************************************************
     *********************************************************************************
     */

    // Index to the first Mutation object not past `start` position:
    uint64 mut_i = var_chrom->get_mut_(start);
    // Current position
    pos = start;

    /*
     If `start` is before the first mutation (resulting in
     `mut_i == var_chrom->mutations.size()`),
     we must pick up any nucleotides before the first mutation.
     */
    if (mut_i == var_chrom->mutations.size()) {
        mut_i = 0;
        for (; pos < var_chrom->mutations[mut_i].new_pos; pos++) {
            cum_wt += nt_rates[var_chrom->ref_chrom->nucleos[pos]];
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
    while (pos <= end && next_mut_i < var_chrom->mutations.size()) {
        while (pos <= end && pos < var_chrom->mutations[next_mut_i].new_pos) {
            char c = var_chrom->get_char_(pos, mut_i);
            cum_wt += nt_rates[c];
            if (cum_wt > u) return;
            ++pos;
        }
        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    while (pos <= end && pos < var_chrom->chrom_size) {
        char c = var_chrom->get_char_(pos, mut_i);
        cum_wt += nt_rates[c];
        if (cum_wt > u) return;
        ++pos;
    }

    for (uint64 i = 0; i < var_chrom->size(); i++) {

        ;

    }

    return;

}


