

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class


#include "mutator_subs.h" // SubMutator
#include "var_classes.h"  // Var* classes
#include "pcg.h"  // runif_01
#include "util.h"  // interrupt_check
#include "alias_sampler.h"  // alias method of sampling


using namespace Rcpp;




int SubMutator::new_rates(const uint64& begin,
                          const uint64& end,
                          std::deque<uint8>& rate_inds,
                          pcg64& eng,
                          Progress& prog_bar) {

    if (!site_var) {
        if (!rate_inds.empty()) rate_inds.clear();
        return 0;
    }

    // (Gammas go from 0 to (n-1), invariants are n.)
    const uint8 n = Q.size();

    const uint64 N = end - begin;
    const uint64 N0 = rate_inds.size();

    uint32 iters = 0;

#ifdef __JACKALOPE_DIAGNOSTICS
    Rcout << std::endl << "rates for " << begin << ' ' << end << ':' << std::endl;
#endif

    if (invariant <= 0) {

        if (N0 > N) rate_inds.resize(N);

        for (uint64 i = 0; i < rate_inds.size(); i++) {
            rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << rate_inds[i] << ' ';
#endif
        }
        while (rate_inds.size() < N) {
            rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << rate_inds.back() << ' ';
#endif
        }

    } else {

        if (N0 > N) rate_inds.resize(N);

        for (uint64 i = 0; i < rate_inds.size(); i++) {
            if (runif_01(eng) > invariant) {
                rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            } else rate_inds[i] = n;
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << rate_inds[i] << ' ';
#endif
        }
        while (rate_inds.size() < N) {
            if (runif_01(eng) > invariant) {
                rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
            } else rate_inds.push_back(n);
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << rate_inds.back() << ' ';
#endif
        }

    }

#ifdef __JACKALOPE_DIAGNOSTICS
    Rcout << std::endl;
#endif

    return 0;
}





inline void SubMutator::adjust_mats(const double& b_len) {

    // UNREST model
    if (U.size() == 0) {
        for (uint32 i = 0; i < Q.size(); i++) {
            // Adjust P(t) matrix using repeated matrix squaring
            Pt_calc(Q[i], 30, b_len, Pt[i]);
            // Now adjust the alias samplers:
            std::vector<AliasSampler>& samp(samplers[i]);
#ifdef __JACKALOPE_DEBUG
            if (samp.size() != 4) stop("SubMutator::adjust_mats-> samp.size() != 4");
#endif
            for (uint32 j = 0; j < 4; j++) {
                samp[j] = AliasSampler(Pt[i].row(j));
            }
        }
    } else {
#ifdef __JACKALOPE_DEBUG
        if (U.size() != Q.size()) stop("SubMutator::adjust_mats-> U.size() != Q.size()");
        if (Ui.size() != Q.size()) stop("SubMutator::adjust_mats-> Ui.size() != Q.size()");
        if (L.size() != Q.size()) stop("SubMutator::adjust_mats-> L.size() != Q.size()");
#endif
        // All other models
        for (uint32 i = 0; i < Q.size(); i++) {
            // Adjust P(t) matrix using eigenvalues and eigenvectors in U, Ui, and L
            Pt_calc(U[i], Ui[i], L[i], b_len, Pt[i]);
            // Now adjust the alias samplers:
            std::vector<AliasSampler>& samp(samplers[i]);
#ifdef __JACKALOPE_DEBUG
            if (samp.size() != 4) stop("SubMutator::adjust_mats-> samp.size() != 4");
#endif
            for (uint32 j = 0; j < 4; j++) {
                samp[j] = AliasSampler(Pt[i].row(j));
            }
        }
    }

    return;

}



//' Add substitutions within a range (pos to (end-1)) before any mutations have occurred.
//'
//' @noRd
//'
inline int SubMutator::subs_before_muts(uint64& pos,
                                        const uint64& begin,
                                        const uint64& end,
                                        const uint8& max_gamma,
                                        const std::string& bases,
                                        const std::deque<uint8>& rate_inds,
                                        VarChrom& var_chrom,
                                        pcg64& eng,
                                        Progress& prog_bar) {

#ifdef __JACKALOPE_DEBUG
    if (rate_inds.empty() && max_gamma > 1) {
        stop("rate_inds shouldn't be empty when max_gamma > 1");
    }
    if (rate_inds.empty() && invariant > 0) {
        stop("rate_inds shouldn't be empty when invariant > 0");
    }
#endif

    uint32 iters = 0;

    if (site_var) {

        for (; pos < end; pos++) {

            const uint8& rate_i(rate_inds[(pos-begin)]);
            if (rate_i > max_gamma) continue; // this is an invariant region

            const uint8& c_i(char_map[var_chrom.ref_chrom->nucleos[pos]]);
            if (c_i > 3) continue; // only changing T, C, A, or G
            AliasSampler& samp(samplers[rate_i][c_i]);
            uint8 nt_i = samp.sample(eng);
            if (nt_i != c_i) var_chrom.add_substitution(bases[nt_i], pos);

            if (interrupt_check(iters, prog_bar)) return -1;

        }


    } else {

        for (; pos < end; pos++) {

            const uint8& c_i(char_map[var_chrom.ref_chrom->nucleos[pos]]);
            if (c_i > 3) continue; // only changing T, C, A, or G
            AliasSampler& samp(samplers.front()[c_i]);
            uint8 nt_i = samp.sample(eng);
            if (nt_i != c_i) var_chrom.add_substitution(bases[nt_i], pos);

            if (interrupt_check(iters, prog_bar)) return -1;

        }

    }


    return 0;

}

//' Add substitutions within a range (pos to (end-1)) after mutations have occurred.
//'
//' @noRd
//'
inline int SubMutator::subs_after_muts(uint64& pos,
                                       const uint64& begin,
                                       const uint64& end1,
                                       const uint64& end2,
                                       const uint64& mut_i,
                                       const uint8& max_gamma,
                                       const std::string& bases,
                                       const std::deque<uint8>& rate_inds,
                                       VarChrom& var_chrom,
                                       pcg64& eng,
                                       Progress& prog_bar) {

    uint64 end = std::min(end1, end2);

    uint32 iters = 0;

    if (site_var) {

        while (pos < end) {

            const uint8& rate_i(rate_inds[(pos-begin)]);
            if (rate_i > max_gamma) {
                pos++;
                continue; // this is an invariant region
            }

            const uint8& c_i(char_map[var_chrom.get_char_(pos, mut_i)]);
            if (c_i > 3) {
                pos++;
                continue; // only changing T, C, A, or G
            }
            AliasSampler& samp(samplers[rate_i][c_i]);
            uint8 nt_i = samp.sample(eng);
            if (nt_i != c_i) var_chrom.add_substitution(bases[nt_i], pos);

            ++pos;

            if (interrupt_check(iters, prog_bar)) return -1;
        }


    } else {

        while (pos < end) {

            const uint8& c_i(char_map[var_chrom.get_char_(pos, mut_i)]);
            if (c_i > 3) {
                pos++;
                continue; // only changing T, C, A, or G
            }
            AliasSampler& samp(samplers.front()[c_i]);
            uint8 nt_i = samp.sample(eng);
            if (nt_i != c_i) var_chrom.add_substitution(bases[nt_i], pos);

            ++pos;

            if (interrupt_check(iters, prog_bar)) return -1;
        }


    }


    return 0;

}






//' Add substitutions for a whole chromosome or just part of one.
//'
//' Here, `end` is NOT inclusive, so can be == var_chrom.size()
//'
//' @noRd
//'
int SubMutator::add_subs(const double& b_len,
                         const uint64& begin,
                         const uint64& end,
                         const std::deque<uint8>& rate_inds,
                         VarChrom& var_chrom,
                         pcg64& eng,
                         Progress& prog_bar) {

#ifdef __JACKALOPE_DEBUG
    if (b_len < 0) {
        Rcout << std::endl << b_len << std::endl;
        stop("b_len < 0 in add_subs");
    }
    if (begin >= var_chrom.size()) {
        Rcout << std::endl << begin << ' ' << var_chrom.size() << std::endl;
        stop("begin >= var_chrom.size() in add_subs");
    }
    if (end > var_chrom.size()) {
        Rcout << std::endl << end << ' ' << var_chrom.size() << std::endl;
        stop("end > var_chrom.size() in add_subs");
    }
#endif


    if ((b_len == 0) || (end == begin)) return 0;

    if (prog_bar.is_aborted() || prog_bar.check_abort()) return -1;

    adjust_mats(b_len);

    uint8 max_gamma = Q.size() - 1; // any rate_inds above this means an invariant region
    std::string bases = "TCAG";

    // To make code less clunky:
    std::deque<Mutation>& mutations(var_chrom.mutations);

    int status;

    /*
     If there are no mutations or if `end-1` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     */
    if (mutations.empty() || ((end-1) < mutations.front().new_pos)) {

        uint64 pos = begin;
        status = subs_before_muts(pos, begin, end, max_gamma, bases, rate_inds, var_chrom,
                                  eng, prog_bar);
        return status;

    }


    // Index to the first Mutation object not past `begin` position:
    uint64 mut_i = var_chrom.get_mut_(begin);
    // Current position
    uint64 pos = begin;

    /*
     If `begin` is before the first mutation (resulting in `mut_i == mutations.size()`),
     we must process any nucleotides before the first mutation.
     */
    if (mut_i == mutations.size()) {

        mut_i = 0;
        status = subs_before_muts(pos, begin, mutations[mut_i].new_pos, max_gamma, bases,
                                  rate_inds, var_chrom, eng, prog_bar);

        if (status < 0) return status;

    }


    /*
     Now, for each subsequent mutation except the last, process all nucleotides
     at or after its position but before the next one.
     */
    uint64 next_mut_i = mut_i + 1;
    while (pos < end && next_mut_i < mutations.size()) {

        status = subs_after_muts(pos, begin, end, mutations[next_mut_i].new_pos, mut_i,
                                 max_gamma, bases, rate_inds, var_chrom, eng, prog_bar);

        if (status < 0) return status;

        ++mut_i;
        ++next_mut_i;
    }

    // Now taking care of nucleotides after the last Mutation
    status = subs_after_muts(pos, begin, end, var_chrom.chrom_size, mut_i,
                             max_gamma, bases, rate_inds, var_chrom, eng, prog_bar);

    return status;

}








// Adjust rate_inds for deletions:
void SubMutator::deletion_adjust(const uint64& size,
                                 uint64 pos,
                                 const uint64& begin,
                                 std::deque<uint8>& rate_inds) {

    if (!site_var) return;

    // Because rate_inds is from `begin` to `end` only
    pos -= begin;

    rate_inds.erase(rate_inds.begin() + pos,
                    rate_inds.begin() + (pos + size));

    return;

}


// Adjust rate_inds for insertions:
void SubMutator::insertion_adjust(const uint64& size,
                                  uint64 pos,
                                  const uint64& begin,
                                  std::deque<uint8>& rate_inds,
                                  pcg64& eng) {

    if (!site_var) return;

    /*
     Because `deque::insert` will insert items before `pos`, and we want it after
     the original `pos`:
     */
    pos++;
    // Because rate_inds is from `begin` to `end` only
    pos -= begin;

    // (Gammas go from 0 to (n-1), invariants are n.)
    const uint8 n = Q.size();

    if (invariant <= 0) {

        for (uint64 i = 0; i < size; i++) {
            rate_inds.insert(rate_inds.begin() + pos,
                             static_cast<uint8>(runif_01(eng) * n));
        }

    } else {

        for (uint64 i = 0; i < size; i++) {
            if (runif_01(eng) > invariant) {
                rate_inds.insert(rate_inds.begin() + pos,
                                 static_cast<uint8>(runif_01(eng) * n));
            } else rate_inds.insert(rate_inds.begin() + pos, n);
        }

    }



    return;
}




