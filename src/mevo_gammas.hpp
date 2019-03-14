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

#include "gemino_types.hpp"       // integer types
#include "pcg.hpp"                // pcg seeding



using namespace Rcpp;



/*
 Stores info for a single Gamma region.
 */

struct GammaRegion {

    double gamma;
    uint32 start;
    uint32 end;

    GammaRegion() {}
    GammaRegion(const double& gamma_, const uint32& start_, const uint32& end_)
        : gamma(gamma_), start(start_), end(end_) {}
    GammaRegion(const GammaRegion& other)
        : gamma(other.gamma), start(other.start), end(other.end) {}
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
    void deletion_adjust(const uint32& ind, std::vector<uint32>& erase_inds,
                         const uint32& del_start, const uint32& del_end,
                         const sint32& del_size);

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

public:

    std::vector<GammaRegion> regions;
    double seq_size;

    SequenceGammas() : regions(), seq_size() {}

    SequenceGammas(const SequenceGammas& other)
        : regions(other.regions), seq_size(other.seq_size) {}

    // Assignment operator
    SequenceGammas& operator=(const SequenceGammas& other) {
        regions = other.regions;
        seq_size = other.seq_size;
        return *this;
    }

    SequenceGammas(arma::mat gamma_mat) {
        if (gamma_mat.n_cols != 2) stop("input Gamma matrix must have 2 columns, "
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
        for (uint32 i = 0; i < gamma_mat.n_rows; i++) {
            // Below, I'm subtracting 1 to go back to 0-based indexing
            regions[i].end = static_cast<uint32>(gamma_mat(i,0)) - 1;
            if (i > 0) {
                regions[i].start = static_cast<uint32>(gamma_mat(i-1,0));
            } else regions[i].start = 0;
            regions[i].gamma = gamma_mat(i,1);
        }
        seq_size = gamma_mat(gamma_mat.n_rows-1, 0);
    }


    /*
     Get Gamma value based on the position on the chromosome.
     */
    inline double operator[](const uint32& pos) const {
        uint32 idx = get_idx(pos);
        return regions[idx].gamma;
    }

    /*
     The same thing as above, except for across a range.
     */
    inline std::vector<double> operator()(const uint32& start, const uint32& end) const {

        std::vector<double> out(end - start + 1);

        uint32 idx = get_idx(start);

        for (uint32 i = start; i <= end;) {
            double gamma = regions[idx].gamma;
            uint32 length = regions[idx].end - i + 1;
            if (i + length - 1 > end) length = end - i + 1;
            std::fill(out.begin() + i - start, out.begin() + i - start + length, gamma);
            idx++;
            i += length;
        }

        return out;
    }

    /*
     Based on a sequence position, return an index to the Gamma region it's inside.
     */
    inline uint32 get_idx(const uint32& pos) const {
        uint32 idx = pos * (static_cast<double>(regions.size()) / seq_size);
        if (idx >= regions.size()) idx = regions.size() - 1;
        while (regions[idx].end < pos) idx++;
        while (regions[idx].start > pos) idx--;
        return idx;
    }


    void update(const uint32& pos, const sint32& size_change);

};



#endif
