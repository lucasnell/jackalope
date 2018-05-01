#ifndef __GEMINO_GAMMAS_H
#define __GEMINO_GAMMAS_H


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution

#include "gemino_types.h" // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "molecular_evolution.h"  //



using namespace Rcpp;




/*
For a sequence, stores all "Gamma regions": regions with the same Gamma modifier to
their overall mutation rate.
Square brackets return a rate of choosing this region: the Gamma value multiplied by its
size. The size was included to account for regions getting larger or smaller due
to indels.
*/


struct GammaRegion {
    double gamma;
    uint start;
    uint end;
    GammaRegion() {}
    GammaRegion(const double& gamma_, const uint& start_, const uint& end_)
        : gamma(gamma_), start(start_), end(end_) {}
};


class SequenceGammas {
private:

    std::vector<GammaRegion> regions;
    double seq_size;

    inline uint get_idx(const uint& new_pos) const {
        if (new_pos > regions.back().end) stop("new_pos too large");
        if (new_pos < regions.front().start) stop("new_pos too small");
        uint idx = new_pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < new_pos) idx++;
        while (regions[idx].start > new_pos) idx--;
        // uint idx = 0;
        // while (regions[idx].end < new_pos) idx++;
        return idx;
    }

    /*
    Inner function to update sizes after a given index.
    */
    inline void update_ranges_after_idx_(uint& idx, const sint& size_change) {
        while (idx < regions.size()) {
            regions[idx].end += size_change;
            regions[idx].start += size_change;
            idx++;
        }
        return;
    }

public:

    arma::mat to_mat() {
        arma::mat out(regions.size(), 3);
        for (uint i = 0; i < regions.size(); i++) {
            out(i,0) = regions[i].start;
            out(i,1) = regions[i].end;
            out(i,2) = regions[i].gamma;
        }
        return out;
    }

    // SequenceGammas(const uint& gamma_size_, const RefSequence& rs,
    //                pcg32& eng, const double& alpha)
    //     : regions(), seq_size(rs.size()) {
    //     // Number of gamma values needed:
    //     uint n_gammas = static_cast<uint>(std::ceil(
    //         static_cast<double>(rs.size()) / static_cast<double>(gamma_size_)));
    //     // Resize gammas and sizes vectors
    //     regions = std::vector<GammaRegion>(n_gammas);
    //
    //     // Fill vectors:
    //     std::gamma_distribution<double> distr(alpha, alpha);
    //     for (uint i = 0, j = 0; i < n_gammas; i++, j+= gamma_size_) {
    //         double gamma_ = distr(eng);
    //         uint end_ = j + gamma_size_ - 1;
    //         if (i == n_gammas - 1) end_ = rs.size() - 1;
    //         regions[i] = GammaRegion(gamma_, j, end_);
    //     }
    // }
    SequenceGammas(const arma::mat& gamma_mat) {
        regions = std::vector<GammaRegion>(gamma_mat.n_rows);
        for (uint i = 0; i < gamma_mat.n_rows; i++) {
            regions[i].start = static_cast<uint>(gamma_mat(i,0));
            regions[i].end = static_cast<uint>(gamma_mat(i,1));
            regions[i].gamma = gamma_mat(i,2);
        }
        seq_size = gamma_mat(gamma_mat.n_rows-1, 1) + 1;
    }

    /*
    Get Gamma value based on the position on the chromosome.
    */
    double get_gamma(const uint& new_pos) const {
        uint idx = get_idx(new_pos);
        return regions[idx].gamma;
    }

    uint size() const noexcept {
        return static_cast<uint>(seq_size);
    }


    /*
    Get relative rate of change compared to other Gamma regions: the Gamma value
    multiplied by the size of the region.
    This index is the position in the `gammas` vector.
    */
    inline double operator[](const uint& idx) const {
        return regions[idx].gamma * (regions[idx].end - regions[idx].start + 1);
    }

    void update_sizes(const uint& new_pos, sint size_change) {

        /*
        -----------
        Do nothing for substitutions
        -----------
        */
        if (size_change == 0) return;

        seq_size += static_cast<double>(size_change);

        // Rcout << "getting index";
        uint idx = get_idx(new_pos);
        // Rcout << " ...done" << std::endl;
        /*
        -----------
        Insertions are also quite simple:
        -----------
        */
        if (size_change > 0) {
            regions[idx].end += size_change;
            idx++;
            update_ranges_after_idx_(idx, size_change); // update all following ranges
            return;
        }

        /*
        -----------
        Deletions... not so much.
        -----------
        */
        const uint& del_start(new_pos);
        uint del_end = new_pos;
        del_end -= (size_change + 1);

        /*
        The simplest (and probably most common) scenario: the deletion is entirely
        within this region but doesn't entirely overlap it.
        */
        if ((del_end <= regions[idx].end && del_start > regions[idx].start) ||
            (del_end < regions[idx].end && del_start >= regions[idx].start)) {
            // Rcout << "simple" << ' ';
            // Rcout << del_end << ' ' << del_start << ", ";
            // Rcout << regions[idx].end << ' ' << regions[idx].start;
            regions[idx].end += size_change;
            idx++;
            update_ranges_after_idx_(idx, size_change); // update all following ranges
            // Rcout << "... done" << std::endl;
            return;
        }

        /*
        Although it shouldn't be common, the scenario where the deletion happens in
        the last region needs to be taken care of now.
        */
        if (idx == regions.size() - 1) {
            // Rcout << "ending" << std::endl;
            // If the deletion overlaps this entire region...
            if (del_end >= regions[idx].end && del_start <= regions[idx].start) {
                regions.erase(regions.begin() + idx);
                return;
            }
            // Else, just adjust the ending position
            regions[idx].end += size_change;
            return;
        }

        /*
        Probably another common scenario: the deletion spans an intersect
        between two regions but doesn't entirely overlap any whole regions
        */
        if (del_end > regions[idx].end && del_end < regions[idx+1].end &&
            del_start > regions[idx].start) {
            // Rcout << "intersection" << std::endl;
            uint overlap = regions[idx].end - del_start + 1;
            regions[idx].end -= overlap;
            regions[idx+1].start -= overlap;
            regions[idx+1].end += size_change;
            idx += 2;
            update_ranges_after_idx_(idx, size_change); // update all following ranges
            return;
        }

        /*
         The last two scenarios should be rare, unless large deletion probabilities
         or sizes are included in simulations.
         Both scenarios involve deletions removing entire region(s).
         */
        // Rcout << "deleting" << std::endl;
        /*
         If deletion doesn't overlap this entire region, take care of its ending
         position, then iterate.
         */
        if (del_start > regions[idx].start) {
            regions[idx].end = del_start - 1;
            idx++;
        }
        // Now erase all regions that this deletion overlaps
        uint start = idx;  // first region to get deleted
        while (del_end >= regions[idx].end) {
            idx++;
            if (idx == regions.size()) break;
        }
        regions.erase(regions.begin() + start, regions.begin() + idx);

        // `start` now points to the point after erased regions

        update_ranges_after_idx_(start, size_change); // update all following ranges

    }



};





#endif
