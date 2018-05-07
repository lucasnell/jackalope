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

struct GammaRegion {

    double gamma;
    uint start;
    uint end;

    GammaRegion() {}
    GammaRegion(const double& gamma_, const uint& start_, const uint& end_)
        : gamma(gamma_), start(start_), end(end_) {}
    // Assignment operator
    GammaRegion& operator=(const GammaRegion& other) {
        gamma = other.gamma;
        start = other.start;
        end = other.end;
        return *this;
    }

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
 Square brackets return the Gamma value for a given sequence position.
 */


class SequenceGammas {

    std::vector<GammaRegion> regions;
    double seq_size;

    /*
     Based on a sequence position, return an index to the Gamma region it's inside.
     */
    inline uint get_idx(const uint& pos) const {
        uint idx = pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < pos) idx++;
        while (regions[idx].start > pos) idx--;
        return idx;
    }

public:

    SequenceGammas(const SequenceGammas& other)
        : regions(other.regions), seq_size(other.seq_size) {}

    SequenceGammas(arma::mat gamma_mat) {
        if (gamma_mat.n_cols != 3) stop("input Gamma matrix must have 2 columns, "
                                            "one for end positions, one for gammas.");
        // Sort from first to last region
        arma::uvec sort_inds = arma::sort_index(gamma_mat.col(0));
        gamma_mat = gamma_mat.rows(sort_inds);
        // Since the input matrix should be 1-based indexing, make sure there are no 0s:
        if (arma::any(gamma_mat.col(0) <= 0)) {
            stop("A value <= 0 was detected in the first column of the Gamma matrix, "
                     "which is where the end points should be. "
                     "Please use 1-based indexing instead of 0-based.");
        }
        // Make sure there are no repeat ends points
        arma::vec diffs = arma::diff(gamma_mat.col(0));
        if (arma::any(diffs == 0)) stop("All Gamma matrix end points must be unique");
        // Now fill in the regions vector
        regions = std::vector<GammaRegion>(gamma_mat.n_rows);
        for (uint i = 0; i < gamma_mat.n_rows; i++) {
            // Below, I'm subtracting 1 to go back to 0-based indexing
            regions[i].end = static_cast<uint>(gamma_mat(i,0)) - 1;
            if (i > 0) {
                regions[i].start = static_cast<uint>(gamma_mat(i-1,0));
            } else regions[i].start = 0;
            regions[i].gamma = gamma_mat(i,1);
        }
        seq_size = gamma_mat(gamma_mat.n_rows-1, 0);
    }
    // Assignment operator
    SequenceGammas& operator=(const SequenceGammas& other) {
        regions = other.regions;
        seq_size = other.seq_size;
        return *this;
    }

    /*
     Get Gamma value based on the position on the chromosome.
     */
    inline double operator[](const uint& pos) const {
        uint idx = get_idx(pos);
        return regions[idx].gamma;
    }

    void update_gamma_regions(const uint& pos, const sint& size_change);

};



arma::mat make_gamma_mat(const uint& seq_size_, const uint& gamma_size_,
                         const double& alpha, pcg32& eng);


#endif
