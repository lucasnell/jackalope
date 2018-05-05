
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
void GammaRegion::deletion_adjust(const uint& ind, std::vector<uint>& erase_inds,
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



void SequenceGammas::update_sizes(const uint& new_pos, sint size_change) {

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
    uint idx = get_idx(new_pos);

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
    const uint& del_start(new_pos);
    uint del_end = new_pos;
    del_end -= (size_change + 1);

    // Iterate through and adjust all regions including and following the deletion:
    std::vector<uint> erase_inds;
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




