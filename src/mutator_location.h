#ifndef __JACKAL_MUTATOR_LOCATION_H
#define __JACKAL_MUTATOR_LOCATION_H


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
#include "site_var.h"  // SequenceGammas class
#include "weighted_reservoir.h"  // weighted_reservoir_* functions
#include "util.h"  // str_stop



using namespace Rcpp;



namespace mut_loc {
    const std::string bases = "TCAG";
}


/*
 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================

 Choosing mutation locations based on overall mutation rates that vary by
 (i) nucleotide and (ii) sequence region

 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================
 */






/*
 Stores info on the overall mutation rates (including Gammas) for each nucleotide.

 This class is used to do weighted reservoir sampling for where a mutation should occur.
 Ultimately this class allows you to use a bracket operator to get a rate for a
 given location in the variant genome.

 For nucleotide rates (not including Gammas), Ns are set to 0 bc we don't want to
 process these.
 Input char objects are cast to uint32 which provide the indices.
 T, C, A, G, and N should never be higher than 84, so will be safe.
 In case someone inputs other characters accidentally, I've set the length to 256,
 which should work for all 8-bit character values.
 The memory overhead should be pretty minimal.
 */

class MutationRates {

public:

    const VarSequence * var_seq;  // pointer to const VarSequence
    std::vector<double> nt_rates;
    SequenceGammas gammas;


    MutationRates() : var_seq(), nt_rates(), gammas() {}

    /*
     Below `q_tcag` is a length-4 vector of rates (q_i from Yang (2006)), for
     T, C, A, and G, respectively
     */
    MutationRates(const VarSequence& vs_, const std::vector<double>& q_tcag,
                  const SequenceGammas& gammas_)
        : var_seq(&vs_), nt_rates(256, 0.0), gammas(gammas_) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 j = mut_loc::bases[i];
            nt_rates[j] = q_tcag[i];
        }
    }
    MutationRates(const std::vector<double>& q_tcag):
        var_seq(), nt_rates(256, 0.0), gammas() {
        for (uint32 i = 0; i < 4; i++) {
            uint32 j = mut_loc::bases[i];
            nt_rates[j] = q_tcag[i];
        }
    }
    MutationRates(const MutationRates& other)
        : var_seq(other.var_seq), nt_rates(other.nt_rates), gammas(other.gammas) {}

    MutationRates& operator=(const MutationRates& other) {
        var_seq = other.var_seq;
        nt_rates = other.nt_rates;
        gammas = other.gammas;
        return *this;
    }


    // To get size of the variant sequence
    // inline uint32 size() const noexcept {
    inline uint32 size() const {
        // If a null pointer, throw error
        if (!var_seq) stop("null pointer when accessing MutationRates");
        return var_seq->size();
    }

    // Return whether the pointer to the VarSequence object is null
    inline bool empty() const noexcept {
        if (!var_seq) return true;
        return false;
    }

    // Using bracket operator to get the overall mutation rate at a location
    inline double operator[](const uint32& pos) const {
        char c = var_seq->get_nt(pos);
        double r = nt_rates[c];
        r *= gammas[pos];
        return r;
    }

    /*
     The same as above, but for a range of positions.
     */
    inline double operator()(const uint32& start, const uint32& end) const {

        std::string seq;
        seq.reserve(end - start + 1);
        uint32 mut_ = var_seq->get_mut_(start);
        var_seq->set_seq_chunk(seq, start, end - start + 1, mut_);

        std::vector<double> gamma_vals = gammas(start, end);
        if (gamma_vals.size() != seq.size()) {
            stop("seq and gamma_vals sizes not matching in MutationRates::().");
        }

        double out = 0;
        for (uint32 i = 0; i < gamma_vals.size(); i++) {
            double r = nt_rates[seq[i]];
            out += r * gamma_vals[i];
        }

        return out;
    }

    /*
     Get the change in mutation rate for a substitution at a location given a
     position and the character it'll change to
     */
    inline double sub_rate_change(const uint32& pos, const char& c) const {
        char c0 = var_seq->get_nt(pos);
        double gamma = gammas[pos];
        double r0 = nt_rates[c0];
        double r1 = nt_rates[c];
        return gamma * (r1 - r0);
    }

    // Return a rate (NO gamma) for an input string
    inline double raw_rate(const std::string& seq) const {
        double out = 0;
        for (const char& c : seq) out += nt_rates[c];
        return out;
    }
    // Overloaded for a character
    inline double raw_rate(const char& seq) const {
        double out = nt_rates[seq];
        return out;
    }

    // Used below to check if gamma region needs to be iterated to the next one.
    inline void check_gamma(const uint32& pos,
                            uint32& gamma_end, uint32& gam_i, double& gamma,
                            const SequenceGammas& gammas) const {
        if (pos > gamma_end) {
            gam_i++;
            gamma = gammas.regions[gam_i].gamma;
            gamma_end = gammas.regions[gam_i].end;
        }
        return;
    }

    // To return the overall rate for an entire sequence:
    double total_rate(uint32 start, uint32 end, const bool& ranged) const {

        double out = 0;

        if (var_seq->size() == 0) return out;

        if (!ranged) {
            start = 0;
            end = var_seq->size() - 1;
        }

        if ((var_seq->size() - 1) != gammas.regions.back().end) {
            stop("gammas and var_seq sizes don't match inside MutationRates");
        }

        /*
         If there are no mutations or if `end` is before the first mutation,
         then we don't need to use the `mutations` field at all.
         (I'm using separate statements to avoid calling `front()` on an empty deque.)
         */
        bool use_mutations = true;
        if (var_seq->mutations.empty()) {
            use_mutations = false;
            if ((var_seq->ref_seq->nucleos.size() - 1) != gammas.regions.back().end) {
                stop("gammas and var_seq ref sizes don't match inside MutationRates");
            }
        } else if (var_seq->mutations.front().new_pos > end) {
            use_mutations = false;
        }
        if (!use_mutations) {

            uint32 i = start, gam_i = gammas.get_idx(start);

            while (i <= end) {
                double gamma = gammas.regions[gam_i].gamma;
                double tmp = 0;
                while (i <= gammas.regions[gam_i].end && i <= end) {
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
        uint32 gam_i = gammas.get_idx(start);

        double gamma = gammas.regions[gam_i].gamma;
        uint32 gamma_end = gammas.regions[gam_i].end;

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
                check_gamma(pos, gamma_end, gam_i, gamma, gammas);
                out += (nt_rates[(*(var_seq->ref_seq))[pos]] * gamma);
            }
            check_gamma(pos, gamma_end, gam_i, gamma, gammas);
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
                check_gamma(pos, gamma_end, gam_i, gamma, gammas);
            }
            ++mut_i;
            ++next_mut_i;
        }

        // Now taking care of nucleotides after the last Mutation
        while (pos <= end &&pos < var_seq->seq_size) {
            char c = var_seq->get_char_(pos, mut_i);
            out += nt_rates[c] * gamma;
            ++pos;
            check_gamma(pos, gamma_end, gam_i, gamma, gammas);
        }

        return out;
    }

};


/*
 This class uses the info above, plus a class and fxn from `weighted_reservoir.h` to
 do weighted reservoir sampling for a single location at which to put a mutation.
 The weights are based on the nucleotide and sequence region.
 */
class LocationSampler {

public:

    ReservoirRates<MutationRates> rates;

    LocationSampler() : rates() {};
    LocationSampler(const MutationRates& mr_, const uint32& chunk = 0)
        : rates(mr_, chunk) {}
    LocationSampler(const LocationSampler& other)
        : rates(other.rates) {}
    LocationSampler& operator=(const LocationSampler& other) {
        rates = other.rates;
        return *this;
    }

    inline uint32 sample(pcg64& eng, const uint32& start, const uint32& end,
                         const bool& ranged) {
        return rates.sample(eng, start, end, ranged);
    }


    inline MutationRates& mr() {
        return rates.res_rates.all_rates;
    }
    inline const MutationRates& mr() const {
        return rates.res_rates.all_rates;
    }

    // Fill pointer
    void fill_ptrs(const VarSequence& vs_) {
        mr().var_seq = &vs_;
        ChunkRateGetter<MutationRates>& chunk_rg(rates.res_rates);
        // Make sure to check on sizes:
        chunk_rg.reset();
        return;
    }

    double substitution_rate_change(const char& c, const uint32& pos) const {
        const MutationRates& mr_(mr());
        return mr_.sub_rate_change(pos, c);
    }

    double insertion_rate_change(const std::string& seq, const uint32& pos) const {
        const MutationRates& mr_(mr());
        double gamma = mr_.gammas[pos];
        double rate = mr_.raw_rate(seq);
        return gamma * rate;
    }

    double deletion_rate_change(const sint32& size_mod, const uint32& start) const {
        uint32 end = start - size_mod - 1;
        const MutationRates& mr_(mr());
        double out = mr_(start, end);
        out *= -1;
        return out;
    }

    inline double total_rate(const uint32& start, const uint32& end,
                             const bool& ranged) const {
        const MutationRates& mr_(mr());
        return mr_.total_rate(start, end, ranged);
    }

    inline void update_gamma_regions(const sint32& size_change, const uint32& pos) {
        MutationRates& mr_(mr());
        mr_.gammas.update(pos, size_change);
        return;
    }


    // Resize chunk size; this method is obviously not available for non-chunked version
    void change_chunk(const uint32& chunk_size) {
        ChunkRateGetter<MutationRates>& crg(rates.res_rates);
        crg.chunk_size = chunk_size;
        crg.inds.resize(chunk_size);
        for (uint32 i = 0; i < chunk_size; i++) crg.inds[i] = i;
        // If not a valid pointer, stop here
        if (!crg.all_rates.var_seq) return;
        // Else, check on sizes:
        if (crg.all_rates.size() > chunk_size) {
            crg.inds.resize(chunk_size);
        } else if (crg.all_rates.size() != crg.inds.size()) {
            crg.inds.resize(crg.all_rates.size());
        }
        // `recheck_size_()` will automatically do the rest of the work from here
        return;
    }
};















#endif
