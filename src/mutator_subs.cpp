

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class


#include "mutator_subs.h" // SubMutator
#include "hap_classes.h"  // Hap* classes
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
        if (!rate_inds.empty()) {
            rate_inds.clear();
        }
        return 0;
    }

    // (Gammas go from 0 to (n-1), invariants are n.)
    const uint8 n = Q.size();

    const uint64 N = end - begin;
    const uint64 N0 = rate_inds.size();

    uint32 iters = 0;

#ifdef __JACKALOPE_DIAGNOSTICS
    Rcout << std::endl << "~~ rates for " << begin << ' ' << end << " = ";
#endif

    if (N0 > N) {
        rate_inds.resize(N);
        clear_memory<std::deque<uint8>>(rate_inds);
    }

    if (invariant <= 0) {

        for (uint64 i = 0; i < rate_inds.size(); i++) {
            rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            // This will be invisible without being converted to unsigned
            Rcout << static_cast<unsigned>(rate_inds[i]) << ' ';
#endif
        }
        while (rate_inds.size() < N) {
            rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << static_cast<unsigned>(rate_inds.back()) << ' ';
#endif
        }

    } else {

        for (uint64 i = 0; i < rate_inds.size(); i++) {
            if (runif_01(eng) > invariant) {
                rate_inds[i] = static_cast<uint8>(runif_01(eng) * n);
            } else rate_inds[i] = n;
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << static_cast<unsigned>(rate_inds[i]) << ' ';
#endif
        }
        while (rate_inds.size() < N) {
            if (runif_01(eng) > invariant) {
                rate_inds.push_back(static_cast<uint8>(runif_01(eng) * n));
            } else rate_inds.push_back(n);
            if (interrupt_check(iters, prog_bar)) return -1;
#ifdef __JACKALOPE_DIAGNOSTICS
            Rcout << static_cast<unsigned>(rate_inds.back()) << ' ';
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
        if (Ui.size() != Q.size()) {
            stop("SubMutator::adjust_mats-> Ui.size() != Q.size()");
        }
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




//' Most of the work for `subs_before_muts`
//'
//' @noRd
//'
inline void SubMutator::subs_before_muts__(const uint64& pos,
                                           uint64& mut_i,
                                           const std::string& bases,
                                           const uint8& rate_i,
                                           HapChrom& hap_chrom,
                                           pcg64& eng) {

    const uint8& c_i(char_map[hap_chrom.ref_chrom->nucleos[pos]]);
    if (c_i > 3) return; // only changing T, C, A, or G
    AliasSampler& samp(samplers[rate_i][c_i]);
    uint8 nt_i = samp.sample(eng);
    if (nt_i != c_i) {
#ifdef __JACKALOPE_DIAGNOSTICS
        // __ <new pos> <rate index> <old nucleotide>-<new nucleotide>
        Rcout << "__ " << pos << ' ' << static_cast<unsigned>(rate_i) << ' ' <<
            bases[c_i] << '-' << bases[nt_i] << std::endl;
#endif
        hap_chrom.mutations.push_front(pos, pos, bases[nt_i]);
        mut_i++;
    }

    return;

}

//' Add substitutions within a range (pos to (end-1)) before any mutations have occurred.
//'
//' @noRd
//'
inline int SubMutator::subs_before_muts(const uint64& begin,
                                        const uint64& end,
                                        uint64& mut_i,
                                        const uint8& max_gamma,
                                        const std::string& bases,
                                        const std::deque<uint8>& rate_inds,
                                        HapChrom& hap_chrom,
                                        pcg64& eng,
                                        Progress& prog_bar,
                                        uint32& iters) {

#ifdef __JACKALOPE_DEBUG
    if (rate_inds.empty() && max_gamma > 1) {
        stop("rate_inds shouldn't be empty when max_gamma > 1");
    }
    if (rate_inds.empty() && invariant > 0) {
        stop("rate_inds shouldn't be empty when invariant > 0");
    }
#endif

    uint64 size = end - begin;
    uint64 pos;

    if (site_var) {

        // Going backwards from `end-1` to `begin` so we can add to front of `mutations`
        for (uint64 i = 1; i <= size; i++) {

            pos = end - i;

            const uint8& rate_i(rate_inds[(pos-begin)]);
            if (rate_i > max_gamma) continue; // this is an invariant region

            subs_before_muts__(pos, mut_i, bases, rate_i, hap_chrom, eng);

            if (interrupt_check(iters, prog_bar)) return -1;

        }


    } else {

        const uint8 rate_i = 0;

        for (uint64 i = 1; i <= size; i++) {

            pos = end - i;

            subs_before_muts__(pos, mut_i, bases, rate_i, hap_chrom, eng);

            if (interrupt_check(iters, prog_bar)) return -1;

        }

    }


    return 0;

}

//' Add substitutions within a range (pos to (end-1)) after mutations have occurred.
//'
//' @noRd
//'
inline void SubMutator::subs_after_muts__(const uint64& pos,
                                          uint64& mut_i,
                                          const std::string& bases,
                                          const uint8& rate_i,
                                          HapChrom& hap_chrom,
                                          pcg64& eng) {


    AllMutations& mutations(hap_chrom.mutations);
    const std::string& reference(hap_chrom.ref_chrom->nucleos);

    const uint8& c_i(char_map[hap_chrom.get_char_(pos, mut_i)]);
    if (c_i > 3) return; // only changing T, C, A, or G

    AliasSampler& samp(samplers[rate_i][c_i]);
    uint8 nt_i = samp.sample(eng);
    const char& nucleo(bases[nt_i]);

    if (nt_i != c_i) {

        sint64 ind = pos - mutations.new_pos[mut_i]; // <-- should always be >= 0

        // If `pos` is within the mutation chromosome:
        if (ind <= hap_chrom.size_modifier(mut_i)) {

#ifdef __JACKALOPE_DIAGNOSTICS
            // __ <new pos> <rate index> <old nucleotide>-<new nucleotide>
            Rcout << "__ " << pos << ' ' << static_cast<unsigned>(rate_i) << ' ' <<
                bases[c_i] << '-' << nucleo << std::endl;
#endif

            /*
             If this new mutation reverts a substitution back to reference state,
             delete the mutation from `mutations`.
             Otherwise, adjust the mutation's sequence.
             When `mut_i == 0`, doing this would make `mut_i` become negative,
             so I just keep the mutation if `mut_i == 0`.
             */
            if ((hap_chrom.size_modifier(mut_i) == 0) &&
                (reference[mutations.old_pos[mut_i]] == nucleo) &&
                mut_i > 0) {
                mutations.erase(mut_i);
                mut_i--;
            } else mutations.nucleos[mut_i][ind] = nucleo;

        } else {
            // If `pos` is in the reference chromosome following the mutation:
            uint64 old_pos_ = ind + (mutations.old_pos[mut_i] -
                hap_chrom.size_modifier(mut_i));
#ifdef __JACKALOPE_DIAGNOSTICS
            // __ <new pos> <rate index> <old nucleotide>-<new nucleotide>
            Rcout << "__ " << pos << ' ' << static_cast<unsigned>(rate_i) << ' ' <<
                bases[c_i] << '-' << nucleo << std::endl;
#endif
            mutations.insert(mut_i + 1, old_pos_, pos, nucleo);
            mut_i++;
        }

    }

    return;

}

//' Add substitutions within a range (pos to (end-1)) after mutations have occurred.
//'
//' @noRd
//'
inline int SubMutator::subs_after_muts(uint64& pos,
                                       const uint64& begin,
                                       const uint64& end1,
                                       const uint64& end2,
                                       uint64& mut_i,
                                       const uint8& max_gamma,
                                       const std::string& bases,
                                       const std::deque<uint8>& rate_inds,
                                       HapChrom& hap_chrom,
                                       pcg64& eng,
                                       Progress& prog_bar,
                                       uint32& iters) {

    uint64 end = std::min(end1, end2);

    if (site_var) {

        while (pos < end) {

            const uint8& rate_i(rate_inds[(pos-begin)]);
            if (rate_i > max_gamma) {
                pos++;
                continue; // this is an invariant region
            }

            subs_after_muts__(pos, mut_i, bases, rate_i, hap_chrom, eng);
            ++pos;
            if (interrupt_check(iters, prog_bar)) return -1;

        }

    } else {

        const uint8 rate_i = 0;

        while (pos < end) {

            subs_after_muts__(pos, mut_i, bases, rate_i, hap_chrom, eng);
            ++pos;
            if (interrupt_check(iters, prog_bar)) return -1;

        }

    }


    return 0;

}







//' Add substitutions for a whole chromosome or just part of one.
//'
//' Here, `end` is NOT inclusive, so can be == hap_chrom.size()
//'
//' @noRd
//'
int SubMutator::add_subs(const double& b_len,
                         const uint64& begin,
                         const uint64& end,
                         const std::deque<uint8>& rate_inds,
                         HapChrom& hap_chrom,
                         pcg64& eng,
                         Progress& prog_bar) {

    if ((b_len == 0) || (end == begin)) return 0;

#ifdef __JACKALOPE_DEBUG
    if (b_len < 0) {
        Rcout << std::endl << b_len << std::endl;
        stop("b_len < 0 in add_subs");
    }
    if (begin >= hap_chrom.size()) {
        Rcout << std::endl << begin << ' ' << hap_chrom.size() << std::endl;
        stop("begin >= hap_chrom.size() in add_subs");
    }
    if (end > hap_chrom.size()) {
        Rcout << std::endl << end << ' ' << hap_chrom.size() << std::endl;
        stop("end > hap_chrom.size() in add_subs");
    }
#endif


    if (prog_bar.is_aborted() || prog_bar.check_abort()) return -1;

    adjust_mats(b_len);

    uint8 max_gamma = Q.size() - 1; // any rate_inds above this means an invariant region
    std::string bases = "TCAG";

    // To make code less clunky:
    AllMutations& mutations(hap_chrom.mutations);

    int status = 0;
    uint32 iters = 0;


    uint64 mut_i = 0;
    uint64 pos = begin;

    /*
     If there are no mutations or if `end-1` is before the first mutation,
     then we don't need to use the `mutations` field at all.
     */
    if (mutations.empty() || ((end-1) < mutations.new_pos.front())) {

        status = subs_before_muts(begin, end, mut_i, max_gamma, bases, rate_inds,
                                  hap_chrom, eng, prog_bar, iters);
        return status;

    }


    /*
     Index to the Mutation object nearest to (without being past) an input position
     on the haplotype chromosome.
     */
    mut_i = hap_chrom.get_mut_(begin);

    /*
     If `begin` is before the first mutation (resulting in `mut_i == mutations.size()`),
     we must process any nucleotides before the first mutation.
     */
    if (mut_i == mutations.size()) {

        mut_i = 0;
        // This is the end for now, but will be `pos` below:
        pos = mutations.new_pos[mut_i];
        status = subs_before_muts(begin, pos, mut_i, max_gamma, bases,
                                  rate_inds, hap_chrom, eng, prog_bar, iters);

        if (status < 0) return status;

    }


    /*
     Now, for each subsequent mutation except the last, process all nucleotides
     at or after its position but before the next one.
     */
    uint64 next_mut_i = mut_i + 1;
    while (pos < end && next_mut_i < mutations.size()) {

        status = subs_after_muts(pos, begin, end, mutations.new_pos[next_mut_i], mut_i,
                                 max_gamma, bases, rate_inds, hap_chrom, eng, prog_bar,
                                 iters);

        if (status < 0) return status;

        ++mut_i;
        next_mut_i = mut_i + 1;
    }

    // Now taking care of nucleotides after the last Mutation
    status = subs_after_muts(pos, begin, end, hap_chrom.chrom_size, mut_i,
                             max_gamma, bases, rate_inds, hap_chrom, eng, prog_bar,
                             iters);

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




