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
                         const sint& del_size);

    inline double size() const {
        return static_cast<double>(end - start + 1);
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

    SequenceGammas(const uint& gamma_size_, const VarSequence& vs,
                   pcg32& eng, const double& alpha)
        : regions(), seq_size(vs.size()) {
        // Number of gamma values needed:
        uint n_gammas = static_cast<uint>(std::ceil(
            static_cast<double>(vs.size()) / static_cast<double>(gamma_size_)));
        // Resize gamma-region vector
        regions = std::vector<GammaRegion>(n_gammas);

        // Fill vector:
        std::gamma_distribution<double> distr(alpha, alpha);
        for (uint i = 0, start_ = 0; i < n_gammas; i++, start_ += gamma_size_) {
            double gamma_ = distr(eng);
            uint end_ = start_ + gamma_size_ - 1;
            if (i == n_gammas - 1) end_ = vs.size() - 1;
            regions[i] = GammaRegion(gamma_, start_, end_);
        }
    }

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
     Notice that this does NOT incorporate the size, so this should
     not be used for sampling.
     */
    double get_gamma(const uint& new_pos) const {
        uint idx = get_idx(new_pos);
        return regions[idx].gamma;
    }

    uint size() const noexcept {
        return static_cast<uint>(seq_size);
    }
    inline double size(const uint& idx) const {
        return regions[idx].size();
    }


    /*
     Get relative rate of change compared to other Gamma regions: the Gamma value
     multiplied by the size of the region.
     This index is the position in the `gammas` vector.
     */
    inline double operator[](const uint& idx) const {
        return regions[idx].gamma * size(idx);
    }

    void update_sizes(const uint& new_pos, sint size_change);

};




#endif
