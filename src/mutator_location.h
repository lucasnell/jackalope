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
struct Region {

    double gamma;
    uint32 start;
    uint32 end;
    double rate;  // overall rate (including gamma) for the sequence in this region

    Region(const double& gamma_, const uint32& start_, const uint32& end_,
                const double& rate_)
        : gamma(gamma_), start(start_), end(end_), rate(rate_) {}
    Region(const Region& other)
        : gamma(other.gamma), start(other.start), end(other.end), rate(other.rate) {}

    /*
     Adjust for a deletion.
     It returns true if this region should be deleted.
     Returns false in all other cases.
     */
    bool del_adjust_bounds(const uint32& del_start,
                           const uint32& del_end);

};



/*
 ========================================================================================
 ========================================================================================

 Binary search tree

 ========================================================================================
 ========================================================================================
 */

// Tips of tree. Mostly just a wrapper around a pointer to a Region.
// This allows the Region to be deleted from memory if its rate == 0.
struct RegionTip {

    std::unique_ptr<Region> region;

    RegionTip() : region(nullptr) {};
    RegionTip(const double& gamma_, const uint32& start_, const uint32& end_,
          const double& rate_)
        : region(new Region(gamma_, start_, end_, rate_)) {}

    inline double rate() {
        if (region) return region->rate;
        return 0;
    }
    // Adjust bounds for an insertion that occurs BEFORE this region:
    void ins_adjust_bounds(const uint32& size) {
        // If not a NULL pointer, adjust the inner Region:
        if (region) {
            region->start += size;
            region->end += size;
        }
        return;
    }
    void del_adjust_bounds(const uint32& del_start,
                           const uint32& del_end) {
        // If not a NULL pointer, adjust the inner Region, and remove the pointer if
        // this deletion covers the whole region:
        if (region) {
            if (region->del_adjust_bounds(del_start, del_end)) {
                region.reset();
            }
        }
        return;
    }

};




// The full tree. This is what outputs a Region position:
class RegionTree {

public:

    std::vector<std::vector<double>> nodes;
    std::vector<RegionTip> tips;
    double total_rate;

    RegionTree(const arma::mat& gamma_mat)
        : total_rate(0) {
        construct_tips(gamma_mat); // also calculates `total_rate`
        construct_nodes();
    }

    /*
     For a sampled rate, return the Region it refers to.
     Note that it changes `u` in place.
     This means that if using CDF method within region, `u` just needs to be divided
     by the region's gamma value to be used.
     */
    Region* search(double& u);

    // This adjusts all node, tip, and rate fields for a substitution:
    void sub_update(const double& d_rate) {
        total_rate += d_rate;
        tips[mut_tip].region->rate += d_rate;
        // Update nodes:
        update_nodes(mut_tip, d_rate);
        // If the region's rate is now zero, remove the object:
        if (tips[mut_tip].region->rate == 0) tips[mut_tip].region.reset();
        return;
    }
    // This adjusts all node, tip, and rate fields for an insertion:
    void ins_update(const double& d_rate, const uint32& size) {
        total_rate += d_rate;
        tips[mut_tip].region->rate += d_rate;
        // Update nodes:
        update_nodes(mut_tip, d_rate);
        // Only update end for the region containing the insertion:
        tips[mut_tip].region->end += size;
        // This method is used for all subsequent regions; updates start and end:
        for (uint32 i = (mut_tip + 1); i < tips.size(); i++) {
            tips[i].ins_adjust_bounds(size);
        }
        return;
    }
    // This adjusts all node, tip, and rate fields for a deletion:
    void del_update(const double& d_rate, const uint32& start, const uint32& end) {
        total_rate += d_rate;
        tips[mut_tip].region->rate += d_rate;
        // Update nodes:
        update_nodes(mut_tip, d_rate);
        // Update region bounds for the one containing the deletion and all
        // subsequent ones:
        for (uint32 i = mut_tip; i < tips.size(); i++) {
            tips[i].del_adjust_bounds(start, end);
        }
        return;
    }

private:

    uint32 mut_tip = 0; // Keeps index to the most recently mutated tip

    void construct_tips(arma::mat gamma_mat);
    void construct_nodes();

    // Given a change in a tip, update all nodes above that for the change
    // It does NOT update the tip itself, NOR does it update the `total_rate` field!
    void update_nodes(const uint32& tip_i, const double& d_rate) {
        if (nodes.size() == 0 || d_rate == 0) return;
        uint32 lvl_i = nodes.size() - 1;
        bool on_left = (tip_i&1) == 0; // using bitwise modulus to see if `tip_i` is even
        uint32 node_i = tip_i>>1; // going up tree means bit-shifting this direction
        // Go up tree and update all necessary nodes:
        while (true) {
            if (on_left) nodes[lvl_i][node_i] += d_rate;
            if (lvl_i == 0) break;
            lvl_i--;
            on_left = (node_i&1) == 0;
            node_i >>= 1;
        }
        return;
    }

};








/*
 ========================================================================================
 ========================================================================================

 Full location sampler class

 ========================================================================================
 ========================================================================================
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
    std::vector<Region> regions;
    double total_rate = 0;
    // For sampling with a starting and ending location:
    uint32 start_pos = 0;
    uint32 end_pos;
    double start_rate = 0;
    double end_rate = 0;
    bool start_end_set = false; // whether pos and rates have been set

    LocationSampler() : var_seq(), regions() {};
    LocationSampler(const VarSequence& vs_,
                    const std::vector<double>& q_tcag,
                    const arma::mat& gamma_mat,
                    const uint32& gamma_size_)
        : var_seq(&vs_), regions(), end_pos(vs_.size()), gamma_size(gamma_size_) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
        construct_gammas(gamma_mat);
    }
    LocationSampler(const std::vector<double>& q_tcag,
                    const uint32& gamma_size_)
        : var_seq(), regions(), gamma_size(gamma_size_) {
        for (uint32 i = 0; i < 4; i++) {
            uint32 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
    }
    LocationSampler(const LocationSampler& other)
        : var_seq(other.var_seq), nt_rates(other.nt_rates),
          regions(other.regions), total_rate(other.total_rate),
          start_pos(other.start_pos), end_pos(other.end_pos),
          start_rate(other.start_rate), end_rate(other.end_rate),
          gamma_size(other.gamma_size) {}
    LocationSampler& operator=(const LocationSampler& other) {
        var_seq = other.var_seq;
        nt_rates = other.nt_rates;
        regions = other.regions;
        total_rate = other.total_rate;
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
        Region& reg(regions[get_gamma_idx(pos)]);
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
        Region& reg(regions[get_gamma_idx(pos)]);
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

    inline void one_gamma_row(const arma::mat& gamma_mat,
                              const uint32& i,
                              uint32& mut_i,
                              std::vector<uint32>& sizes);


    // Inner method that does most of the work for `calc_rate`
    double calc_rate__(uint32 start, uint32 end) const;


    inline void gamma_sample(uint32& pos,
                             double& u,
                             double& cum_wt,
                             const uint32& gam_i) const;

    inline void safe_get_mut(const uint32& pos, uint32& mut_i) const;


    inline double partial_gamma_rate___(const uint32& end,
                                        const Region& reg) const;


    void update_start_end(const uint32& start, const uint32& end);

};















#endif
