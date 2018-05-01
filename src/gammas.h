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


class GammaRegion {
public:
    double gamma;
    uint start;
    uint end;
    GammaRegion() {}
    GammaRegion(const double& gamma_, const uint& start_, const uint& end_)
        : gamma(gamma_), start(start_), end(end_) {}

    /*
     Adjust for a deletion.
     `ind` is the index to the current region in the vector of regions.
     `erase_inds` stores indices for region(s) to be erased if the deletion
     entirely spans one or more region(s).
     Adding to this variable will result in the current region being erased.
     */
    void deletion_adjust(const uint& ind, std::vector<uint>& erase_inds,
                         const uint& del_start, const uint& del_end,
                         const sint& del_size) {

        // Total overlap
        if (del_start <= start & del_end >= end) {
            erase_inds.push_back(ind);
            return;
        }
        // Deletion is totally inside this region but doesn't entirely overlap it
        if (del_start > start & del_end < end) {
            end += del_size;
            return;
        }
        // No overlap and deletion starts after it
        if (del_start > end) return;
        // No overlap and deletion starts before it
        if (del_end < start) {
            start += del_size;
            end += del_size;
            return;
        }
        // Partial overlap at the start
        if (del_end >= start & del_start < start) {
            start = del_start;
            end += del_size;
            return;
        }
        // Partial overlap at the end
        if (del_start <= end & del_end > end) {
            end = del_start - 1;
        }
        return;
    }
};


class SequenceGammas {
private:

    std::vector<GammaRegion> regions;
    double seq_size;

    inline uint get_idx(const uint& new_pos) const {
        uint idx = new_pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < new_pos) idx++;
        while (regions[idx].start > new_pos) idx--;
        return idx;
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

        uint idx = get_idx(new_pos);

        /*
         -----------
         Insertions are also quite simple:
         -----------
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

        std::vector<uint> erase_inds;
        while (idx < regions.size()) {
            regions[idx].deletion_adjust(idx, erase_inds, del_start, del_end,
                                         size_change);
            idx++;
        }

        // If any regions need erasing, their indices will be stored in erase_inds
        if (erase_inds.size() == 1) {
            regions.erase(regions.begin() + erase_inds.front());
        } else if (erase_inds.size() > 1) {
            regions.erase(regions.begin() + erase_inds.front(),
                          regions.begin() + erase_inds.back() + 1);
        }

    }

};




#endif
