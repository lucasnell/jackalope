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
     If the deletion totally overlaps it, this function changes this regions's rate
     to -1, and makes both start and end 0;
     it lastly returns true.
     Otherwise, it returns false.
    */
    bool deletion_adjust(const uint32& del_start,
                         const uint32& del_end,
                         const uint32& del_size);

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
    long double total_rate = 0;
    /*
     Whether to sample within a GammaRegion using rejection sampling.
     The alternative is to make GammaRegions quite small and use CDF method.
     */
    bool rej_sample;
    double max_q;
    // For sampling with a starting and ending location:
    uint32 start_pos = 0;
    uint32 end_pos;
    long double start_rate = 0;
    long double end_rate = 0;
    bool start_end_set = false; // whether pos and rates have been set

    // If gamma_size == 0, then it won't split GammaRegions
    LocationSampler() : var_seq(), regions() {};
    LocationSampler(const VarSequence& vs_,
                    const std::vector<double>& q_tcag,
                    const bool& rej_sample_,
                    const arma::mat& gamma_mat,
                    const uint32& gamma_size_)
        : var_seq(&vs_), regions(),
          rej_sample(rej_sample_),
          max_q(*std::max_element(q_tcag.begin(), q_tcag.end())),
          end_pos(vs_.size()), gamma_size(gamma_size_) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }

        construct_gammas(gamma_mat);
    }
    LocationSampler(const std::vector<double>& q_tcag,
                    const bool& rej_sample_,
                    const uint32& gamma_size_)
        : var_seq(), regions(),
          rej_sample(rej_sample_),
          max_q(*std::max_element(q_tcag.begin(), q_tcag.end())),
          gamma_size(gamma_size_) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
    }
    LocationSampler(const LocationSampler& other)
        : var_seq(other.var_seq), nt_rates(other.nt_rates),
          regions(other.regions), total_rate(other.total_rate),
          rej_sample(other.rej_sample), max_q(other.max_q),
          start_pos(other.start_pos), end_pos(other.end_pos),
          start_rate(other.start_rate), end_rate(other.end_rate),
          gamma_size(other.gamma_size) {}
    LocationSampler& operator=(const LocationSampler& other) {
        var_seq = other.var_seq;
        nt_rates = other.nt_rates;
        regions = other.regions;
        total_rate = other.total_rate;
        rej_sample = other.rej_sample;
        max_q = other.max_q;
        start_pos = other.start_pos;
        end_pos = other.end_pos;
        start_rate = other.start_rate;
        end_rate = other.end_rate;
        gamma_size = other.gamma_size;
        return *this;
    }


    uint32 sample(pcg64& eng, const uint32& start, const uint32& end);
    uint32 sample(pcg64& eng) const;


    // Fill pointer for a new VarSequence
    void new_seq(const VarSequence& vs_, const arma::mat& gamma_mat) {
        var_seq = &vs_;
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
        end_rate += d_rate;
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
        end_rate += d_rate;
        return d_rate;
    }

    double deletion_rate_change(const uint32& del_size, const uint32& start);

    double calc_rate() const;
    double calc_rate(const uint32& start, const uint32& end) const;

    void update_gamma_regions(const sint32& size_change,
                              const uint32& pos);



private:

    uint32 gamma_size;


    void construct_gammas(arma::mat gamma_mat);


    /*
     Based on a sequence position, return an index to the Gamma region it's inside.
     */
    inline uint32 get_gamma_idx(const uint32& pos) const {
        uint32 idx = pos * (static_cast<double>(regions.size()) /
            static_cast<double>(var_seq->size()));
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < pos || regions[idx].rate < 0) idx++;
        while (regions[idx].start > pos || regions[idx].rate < 0) idx--;
        return idx;
    }

    // Used to check if gamma region needs to be iterated to the next one.
    inline void check_gamma(const uint32& pos,
                            uint32& gamma_end,
                            uint32& gam_i,
                            double& gamma) const;

    inline void one_gamma_row(const arma::mat& gamma_mat,
                              const uint32& i,
                              uint32& mut_i,
                              std::vector<uint32>& sizes);


    // Inner method that does most of the work for `calc_rate`
    double calc_rate__(const uint32& start, const uint32& end) const;


    // Sample within one GammaRegion using rejection method:
    inline void rej_region_sample(uint32& pos,
                                  const uint32& start_,
                                  const uint32& end_,
                                  pcg64& eng,
                                  const uint32& gam_i) const;
    // Sample within one GammaRegion using cdf (non-rejection) method:
    inline void cdf_region_sample(uint32& pos,
                                  long double& u,
                                  long double& cum_wt,
                                  const uint32& gam_i) const;

    inline void safe_get_mut(const uint32& pos, uint32& mut_i) const;


    inline long double partial_gamma_rate___(const uint32& end,
                                             const GammaRegion& reg) const;


    void update_start_end(const uint32& start, const uint32& end);

};















#endif
