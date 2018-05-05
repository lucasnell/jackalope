#ifndef __GEMINO_GAMMAS_H
#define __GEMINO_GAMMAS_H


/*
 ****************************
 ****************************

 Classes related to Gamma regions: sequence regions with the same Gamma modifier to
 the overall mutation rate.

 ****************************
 ****************************
 */

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp>   // pcg prng
#include <vector>               // vector class
#include <string>               // string class
#include <random>               // gamma_distribution

#include "gemino_types.h"       // integer types
#include "sequence_classes.h"   // Var* and Ref* classes
#include "pcg.h"                // pcg seeding



using namespace Rcpp;



/*
 Stores info for a single Gamma region.
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



/*
 For a sequence, stores all "Gamma regions": regions with the same Gamma modifier to
 their overall mutation rate.
 Square brackets return a rate of choosing this region: the Gamma value multiplied by its
 size. The size was included to account for regions getting larger or smaller due
 to indels.
 */


class SequenceGammas {
private:

    std::vector<GammaRegion> regions;
    double seq_size;

    /*
     Based on a sequence position, return an index to the Gamma region it's inside.
     */
    inline uint get_idx(const uint& new_pos) const {
        uint idx = new_pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < new_pos) idx++;
        while (regions[idx].start > new_pos) idx--;
        return idx;
    }


public:

    SequenceGammas(const SequenceGammas& other)
        : regions(other.regions), seq_size(other.seq_size) {}

    SequenceGammas(const VarSequence& vs, const uint& gamma_size_,
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

    SequenceGammas(arma::mat gamma_mat) {
        // Sort from first to last region
        arma::uvec sort_inds = arma::sort_index(gamma_mat.col(0));
        gamma_mat = gamma_mat.rows(sort_inds);
        // Check that matrix regions are appropriately labelled
        std::string err = "input Gamma matrix needs to use 0-based, inclusive indexing, ";
        err += "have the starting points in the first column, ";
        err += "and have the ending points in the second column.";
        if (gamma_mat(0,0) != 0) stop(err);
        for (uint i = 1; i < gamma_mat.n_rows; i++) {
            if ((gamma_mat(i,1) - gamma_mat(i-1,1)) != 1) stop(err);
        }
        // Now fill in the regions vector
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
    inline double operator[](const uint& new_pos) const {
        uint idx = get_idx(new_pos);
        return regions[idx].gamma;
    }

    void update_sizes(const uint& new_pos, sint size_change);

};




#endif
