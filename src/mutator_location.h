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
 Stores info for a single Gamma region.
 This struct doesn't do much other than store the basic info.
 */

struct GammaRegion {

    double gamma;
    uint32 start;
    uint32 end;
    long double rate;  // overall rate (including gamma) for the sequence in this region

    GammaRegion() {}
    GammaRegion(const double& gamma_, const uint32& start_, const uint32& end_,
                const long double& rate_)
        : gamma(gamma_), start(start_), end(end_), rate(rate_) {}
    GammaRegion(const GammaRegion& other)
        : gamma(other.gamma), start(other.start), end(other.end), rate(other.rate) {}
    // Assignment operator
    GammaRegion& operator=(const GammaRegion& other) {
        gamma = other.gamma;
        start = other.start;
        end = other.end;
        rate = other.rate;
        return *this;
    }

    /*
    Adjust for a deletion.
    `ind` is the index to the current region in the vector of regions.
    `erase_inds` stores indices for region(s) to be erased if the deletion
    entirely spans one or more region(s).
    Adding to this variable will result in the current region being erased.
    */
    void deletion_adjust(const uint32& ind,
                         std::vector<uint32>& erase_inds,
                         const uint32& del_start,
                         const uint32& del_end,
                         const sint32& del_size);

    inline double size() const {
        return static_cast<double>(end - start + 1);
    }

};





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





/*
 This class uses the info above, plus a class and fxn from `weighted_reservoir.h` to
 do weighted reservoir sampling for a single location at which to put a mutation.
 The weights are based on the nucleotide and sequence region.
 */
class LocationSampler {

public:

    const VarSequence * var_seq;  // pointer to const VarSequence
    std::vector<double> nt_rates = std::vector<double>(256, 0.0);
    std::vector<GammaRegion> regions;
    double seq_size;
    long double total_rate = 0;

    LocationSampler() : var_seq(), regions(), seq_size() {};
    LocationSampler(const VarSequence& vs_,
                    const std::vector<double>& q_tcag,
                    const arma::mat& gamma_mat)
        : var_seq(&vs_), regions(), seq_size(vs_.size()) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
        construct_gammas(gamma_mat);
    }
    LocationSampler(const std::vector<double>& q_tcag)
        : var_seq(), regions(), seq_size() {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
    }
    LocationSampler(const LocationSampler& other)
        : var_seq(other.var_seq), nt_rates(other.nt_rates),
          regions(other.regions), seq_size(other.seq_size),
          total_rate(other.total_rate) {}
    LocationSampler& operator=(const LocationSampler& other) {
        var_seq = other.var_seq;
        nt_rates = other.nt_rates;
        regions = other.regions;
        seq_size = other.seq_size;
        total_rate = other.total_rate;
        return *this;
    }

    // inline uint32 size() const noexcept {
    inline uint32 size() const {
        // If a null pointer, throw error
        if (!var_seq) stop("null pointer when accessing LocationSampler");
        return var_seq->size();
    }

    uint32 sample(pcg64& eng, const uint32& start, const uint32& end);
    uint32 sample(pcg64& eng);


    // Fill pointer for a new VarSequence
    void new_seq(const VarSequence& vs_, const arma::mat& gamma_mat) {
        var_seq = &vs_;
        seq_size = vs_.size();
        construct_gammas(gamma_mat);
        return;
    }

    double substitution_rate_change(const char& c, const uint32& pos) {
        GammaRegion& reg(regions[get_gamma_idx(pos)]);
        char c0 = var_seq->get_nt(pos);
        double gamma = reg.gamma;
        double d_rate = nt_rates[c] - nt_rates[c0];
        d_rate *= gamma;
        reg.rate += d_rate;
        total_rate += d_rate;
        return d_rate;
    }

    double insertion_rate_change(const std::string& seq, const uint32& pos) {
        GammaRegion& reg(regions[get_gamma_idx(pos)]);
        double gamma = reg.gamma;
        double d_rate = 0;
        for (const char& c : seq) d_rate += nt_rates[c];
        d_rate *= gamma;
        reg.rate += d_rate;
        total_rate += d_rate;
        return d_rate;
    }

    double deletion_rate_change(const sint32& size_mod, const uint32& start);

    double calc_rate() const;
    double calc_rate(const uint32& start, const uint32& end) const;

    void update_gamma_regions(const sint32& size_change,
                              const uint32& pos);





private:


    void construct_gammas(arma::mat gamma_mat) {

        total_rate = 0;

        if (!var_seq) stop("Cannot do construct_gammas method when var_seq isn't set.");
        // Sort from first to last region
        arma::uvec sort_inds = arma::sort_index(gamma_mat.col(0));
        gamma_mat = gamma_mat.rows(sort_inds);
        // Now fill in the regions vector
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
                start = static_cast<uint32>(gamma_mat(i-1,0));
            } else start = 0;

            gamma = gamma_mat(i,1);

            std::string seq;
            rate = 0;
            var_seq->set_seq_chunk(seq, start, end - i + 1, mut_i);
            for (const char& c : seq) rate += nt_rates[c];

            regions.push_back(GammaRegion(gamma, start, end, rate));

            total_rate += rate;
        }

        return;
    }


    /*
     Based on a sequence position, return an index to the Gamma region it's inside.
     */
    inline uint32 get_gamma_idx(const uint32& pos) const {
        uint32 idx = pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < pos) idx++;
        while (regions[idx].start > pos) idx--;
        return idx;
    }

    // Used to check if gamma region needs to be iterated to the next one.
    inline void check_gamma(const uint32& pos,
                            uint32& gamma_end,
                            uint32& gam_i,
                            double& gamma) const {
        if (pos > gamma_end) {
            gam_i++;
            gamma = regions[gam_i].gamma;
            gamma_end = regions[gam_i].end;
        }
        return;
    }


    // Inner method that does most of the work for `calc_rate`
    double calc_rate__(uint32 start, uint32 end) const;

};















#endif
