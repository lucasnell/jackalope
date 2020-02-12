
/*
 Combining samplers for substitutions and indels into a mutation sampler for
 evolving chromosomes along trees.
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "mutator.h"   // TreeMutator
#include "mutator_subs.h"   // SubMutator
#include "mutator_indels.h" // IndelMutator
#include "io.h"  // FileUncomp
#include "util.h"  // str_stop



using namespace Rcpp;



/*
 Add mutations for a branch within a range.
 It also updates `end` for indels that occur in the range.
 (`end == begin` when chromosome region is of size zero bc `end` is non-inclusive)
 */
int TreeMutator::mutate(const double& b_len,
                        HapChrom& hap_chrom,
                        pcg64& eng,
                        Progress& prog_bar,
                        const uint64& begin,
                        uint64& end,
                        std::deque<uint8>& rate_inds)  {

#ifdef __JACKALOPE_DEBUG
    if (end < begin) stop("end < begin in TreeMutator.mutate");
    if (end == begin) stop("end == begin in TreeMutator.mutate");
#endif

    int status;

    status = indels.add_indels(b_len, begin, end, rate_inds, subs, hap_chrom,
                               eng, prog_bar);
    if (status < 0) return status;

    status = subs.add_subs(b_len, begin, end, rate_inds, hap_chrom, eng, prog_bar);

    return status;
}

