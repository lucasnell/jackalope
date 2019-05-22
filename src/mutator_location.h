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
    uint64 start;
    uint64 end;
    double rate;  // overall rate (including gamma) for the sequence in this region
    bool deleted;

    Region(const double& gamma_, const uint64& start_, const uint64& end_,
                const double& rate_)
        : gamma(gamma_), start(start_), end(end_), rate(rate_), deleted(rate_ <= 0) {}
    Region(const Region& other)
        : gamma(other.gamma), start(other.start), end(other.end), rate(other.rate),
          deleted(other.deleted) {}

    /*
     Adjust for a deletion.
     */
    void del_adjust_bounds(const uint64& del_start,
                           const uint64& del_end);

};



/*
 ========================================================================================
 ========================================================================================

 Binary search tree

 ========================================================================================
 ========================================================================================
 */



// The full tree. This is what outputs a Region position:
class RegionTree {

public:

    std::vector<std::vector<double>> nodes;
    std::vector<Region> tips;
    double total_rate;

    RegionTree() {};
    RegionTree(const arma::mat& gamma_mat,
               const uint64& region_size,
               const VarSequence* var_seq,
               const std::vector<double>& nt_rates)
        : total_rate(0) {
        // below also calculates `total_rate`
        construct_tips(gamma_mat, region_size, var_seq, nt_rates);
        construct_nodes();
    }
    RegionTree(const RegionTree& other)
        : nodes(other.nodes), tips(other.tips), total_rate(other.total_rate) {}
    RegionTree& operator=(const RegionTree& other) {
        nodes = other.nodes;
        tips = other.tips;
        total_rate = other.total_rate;
        return *this;
    }

    /*
     For a sampled number in range (0, `total_rate`) or
     (`bounds.start_rate`,`bounds.end_rate`), return the Region it refers to.
     This function adjusts `u` in place so that it refers to a location inside
     the region returned.
     This means that if using CDF method within region, `u` just needs to be divided
     by the region's gamma value to be used.
     */
    const Region* search(double& u) const;

    const Region* current() const {
        return &tips[mut_tip_];
    }

    // Get index to mutated tip (but don't allow outside classes to modify it)
    const uint64& mut_tip() const { return mut_tip_; }

    // This adjusts all node, tip, and rate fields for a substitution:
    void sub_update(const double& d_rate) {
        total_rate += d_rate;
        tips[mut_tip_].rate += d_rate;
        // Update nodes:
        update_nodes(mut_tip_, d_rate);
        // If the region's rate is now zero, remove the object:
        if (tips[mut_tip_].rate == 0) tips[mut_tip_].deleted = true;
        return;
    }
    // This adjusts all node, tip, and rate fields for an insertion:
    void ins_update(const double& d_rate, const uint64& size) {
        total_rate += d_rate;
        tips[mut_tip_].rate += d_rate;
        // Update nodes:
        update_nodes(mut_tip_, d_rate);
        // Only update end for the region containing the insertion:
        tips[mut_tip_].end += size;
        // This method is used for all subsequent regions; updates start and end:
        for (uint64 i = (mut_tip_ + 1); i < tips.size(); i++) {
            tips[i].start += size;
            tips[i].end += size;
        }
        return;
    }
    // This adjusts all node, tip, and rate fields for a deletion:
    void del_update(const double& d_rate,
                    const uint64& start,
                    const uint64& end,
                    std::deque<double>& del_rate_changes) {
        total_rate += d_rate;
        /*
         We don't want to do the following 2 lines bc deletion rate changes are saved
         inside `del_rate_changes`. This is done bc `d_rate` can span multiple regions.
         */
        // tips[mut_tip_].rate += d_rate;
        // update_nodes(mut_tip_, d_rate);

        // Update nodes:
        // Update region bounds for the one containing the deletion and all
        // subsequent ones:
        for (uint64 i = mut_tip_; i < tips.size(); i++) {
            tips[i].del_adjust_bounds(start, end);
            // If `del_rate_changes` isn't empty, then this region was affected:
            if (!del_rate_changes.empty()) {
                double d_rate_ = del_rate_changes.front();
                tips[i].rate += d_rate_;
                update_nodes(i, d_rate_);
                del_rate_changes.pop_front();
            }
        }
        return;
    }

private:

    mutable uint64 mut_tip_ = 0; // Keeps index to the most recently mutated tip

    inline void construct_tips_one_row(const arma::mat& gamma_mat,
                                       const uint64& region_size,
                                       const VarSequence* var_seq,
                                       const std::vector<double>& nt_rates,
                                       const uint64& i,
                                       uint64& mut_i,
                                       std::vector<uint64>& sizes);
    void construct_tips(arma::mat gamma_mat,
                        const uint64& region_size,
                        const VarSequence* var_seq,
                        const std::vector<double>& nt_rates);
    void construct_nodes();

    // Given a change in a tip, update all nodes above that for the change
    // It does NOT update the tip itself, NOR does it update the `total_rate` field!
    void update_nodes(const uint64& tip_i, const double& d_rate) {
        if (nodes.size() == 0 || d_rate == 0) return;
        uint64 lvl_i = nodes.size() - 1;
        bool on_left = (tip_i&1) == 0; // using bitwise modulus to see if `tip_i` is even
        uint64 node_i = tip_i>>1; // going up tree means bit-shifting this direction
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

 For setting bounds to sampling

 ========================================================================================
 ========================================================================================
 */


struct LocationBounds {

    uint64 start_pos;
    uint64 end_pos;
    double start_rate;
    double end_rate;
    bool start_end_set; // whether pos and rates have been set

    LocationBounds()
        : start_pos(0), end_pos(0),
          start_rate(0), end_rate(0),
          start_end_set(false) {};

    LocationBounds(const VarSequence& vs_)
        : start_pos(0), end_pos(vs_.size()),
          start_rate(0), end_rate(0),
          start_end_set(false) {};

    LocationBounds(const LocationBounds& other)
        : start_pos(other.start_pos), end_pos(other.end_pos),
          start_rate(other.start_rate), end_rate(other.end_rate),
          start_end_set(other.start_end_set) {};

    LocationBounds& operator=(const LocationBounds& other) {
        start_pos = other.start_pos;
        end_pos = other.end_pos;
        start_rate = other.start_rate;
        end_rate = other.end_rate;
        start_end_set = other.start_end_set;
        return *this;
    }

    void update(const double& d_rate, const sint64& size_change) {
        end_rate += d_rate;
        end_pos += size_change;
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



class LocationSampler {

public:

    const VarSequence * var_seq;  // pointer to const VarSequence
    std::vector<double> nt_rates = std::vector<double>(256, 0.0);
    // Binary search tree for sampling regions:
    RegionTree regions;
    // For sampling with a starting and ending location:
    LocationBounds bounds;

    LocationSampler() : var_seq(), regions(), bounds(), region_size(0) {};
    LocationSampler(const std::vector<double>& q_tcag,
                    const uint64& region_size_)
        : var_seq(), regions(), bounds(), region_size(region_size_) {
        for (uint64 i = 0; i < 4; i++) {
            uint64 bi = mut_loc::bases[i];
            nt_rates[bi] = q_tcag[i];
        }
    }
    LocationSampler(const LocationSampler& other)
        : var_seq(other.var_seq), nt_rates(other.nt_rates), regions(other.regions),
          bounds(other.bounds), region_size(other.region_size) {}
    LocationSampler& operator=(const LocationSampler& other) {
        var_seq = other.var_seq;
        nt_rates = other.nt_rates;
        regions = other.regions;
        bounds = other.bounds;
        region_size = other.region_size;
        return *this;
    }


    uint64 sample(pcg64& eng, const uint64& start, const uint64& end);
    uint64 sample(pcg64& eng) const;


    // Fill pointer for a new VarSequence
    void new_seq(const VarSequence& vs_, const arma::mat& gamma_mat) {
        var_seq = &vs_;
        // Reset search tree for regions:
        regions = RegionTree(gamma_mat, region_size, var_seq, nt_rates);
        // Reset bounds:
        bounds = LocationBounds(vs_);
        return;
    }

    double total_rate() const { return regions.total_rate; }


    /*
     These figure out the change in rates for each type of mutation, but they do
     NOT change anything.
     */

    double substitution_rate_change(const char& c, const uint64& pos) const {
        // Pointer to region where mutation occurred (saved by `regions` object)
        const Region* reg = regions.current();
        char c0 = var_seq->get_nt(pos);
        double d_rate = nt_rates[c] - nt_rates[c0];
        d_rate *= reg->gamma;
        return d_rate;
    }

    double insertion_rate_change(const std::string& seq, const uint64& pos) const {
        const Region* reg = regions.current();
        double d_rate = 0;
        for (const char& c : seq) d_rate += nt_rates[c];
        d_rate *= reg->gamma;
        return d_rate;
    }

    double deletion_rate_change(const uint64& del_size, const uint64& start) const;



    /*
     These DO update and change things:
     */
    void update(const double& d_rate, const sint64& size_change, const uint64& pos);

    void new_bounds(const uint64& start, const uint64& end);



private:

    uint64 region_size;
    // to store info on deletion rate changes since they affect multiple regions:
    mutable std::deque<double> del_rate_changes;


    // Sample within a region using CDF method:
    inline void cdf_region_sample(uint64& pos, double& u, const Region* reg) const;

    inline void safe_get_mut(const uint64& pos, uint64& mut_i) const;


    inline double partial_gamma_rate___(const uint64& end,
                                        const Region& reg) const;

};















#endif
