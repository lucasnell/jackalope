
/*
 ****************************
 ****************************

 Methods for classes related to Gamma regions: sequence regions with the same Gamma
 modifier to the overall mutation rate.

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
#include "mevo_gammas.h"        // Gamma* classes



using namespace Rcpp;




/*
 Adjust for a deletion.
 `ind` is the index to the current region in the vector of regions.
 `erase_inds` stores indices for region(s) to be erased if the deletion
 entirely spans one or more region(s).
 Adding to this variable will result in the current region being erased.
 */
void GammaRegion::deletion_adjust(const uint32& ind, std::vector<uint32>& erase_inds,
                                  const uint32& del_start, const uint32& del_end,
                                  const sint32& del_size) {

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



void SequenceGammas::update(const uint32& pos, const sint32& size_change) {

    /*
     -----------
     Substitutions
     -----------
     */
    if (size_change == 0) return;


    /*
     -----------
     InDels
     -----------
     */

    seq_size += static_cast<double>(size_change);
    uint32 idx = get_idx(pos);


    /*
     Insertions
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
     Deletions
     */
    const uint32& del_start(pos);
    uint32 del_end = pos;
    del_end -= (size_change + 1);

    // Iterate through and adjust all regions including and following the deletion:
    std::vector<uint32> erase_inds;
    while (idx < regions.size()) {
        regions[idx].deletion_adjust(idx, erase_inds, del_start, del_end,
                                     size_change);
        idx++;
    }

    /*
     If any regions need erasing, their indices will be stored in erase_inds.
     They will be consecutive indices, so we only need to access the front and back.
     */
    if (erase_inds.size() == 1) {
        regions.erase(regions.begin() + erase_inds.front());
    } else if (erase_inds.size() > 1) {
        regions.erase(regions.begin() + erase_inds.front(),
                      regions.begin() + erase_inds.back() + 1);
    }

    return;
}




//' Make matrix of Gamma-region end points and Gamma values.
//'
//' @param seq_size_ Length of the focal sequence.
//' @param gamma_size_ Size of each Gamma region.
//' @param alpha The alpha (shape) parameter for the Gamma distribution from which
//'     Gamma values will be derived.
//' @param eng A random number generator.
//'
//'
//' @noRd
//'
arma::mat make_gamma_mat(const uint32& seq_size_, const uint32& gamma_size_,
                         const double& alpha, pcg32& eng) {

    // If gamma_size_ is set to zero, then we'll assume everything's the same
    if (gamma_size_ == 0) {
        arma::mat out(1, 2, arma::fill::zeros);
        out(0,0) = seq_size_;
        out(0,1) = 1;
        return out;
    }

    // Number of gamma values needed:
    uint32 n_gammas = static_cast<uint32>(std::ceil(
        static_cast<double>(seq_size_) / static_cast<double>(gamma_size_)));

    // Initialize output matrix
    arma::mat out(n_gammas, 2, arma::fill::zeros);

    // Initialize Gamma distribution
    std::gamma_distribution<double> distr(alpha, alpha);

    /*
     Fill matrix.
     Note that I'm not subtracting 1 bc the SequenceGammas constructor assumes
     1-based indexing.
     I'm doing it this way to make it more straightforward if someone wants to pass
     their own matrix directly from R (since R obviously uses 1-based).
     */
    for (uint32 i = 0, start_ = 0; i < n_gammas; i++, start_ += gamma_size_) {
        double gamma_ = distr(eng);
        uint32 end_ = start_ + gamma_size_;
        if (i == n_gammas - 1) end_ = seq_size_;
        out(i,0) = end_;
        out(i,1) = gamma_;
    }

    return out;
}


