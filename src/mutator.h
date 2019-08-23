#ifndef __JACKALOPE_MUTATOR_H
#define __JACKALOPE_MUTATOR_H


/*
 Combining samplers for substitutions and indels into a mutation sampler for
 a single chromosome.
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "jackalope_types.h" // integer types
#include "var_classes.h"  // Var* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "mutator_subs.h"   // SubMutator
#include "mutator_indels.h" // IndelMutator
#include "io.h"  // FileUncomp
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 MutationSampler combines objects for sampling substitutions and indels.
 */
class MutationSampler {

public:

    // VarChrom object pointer to be manipulated
    VarChrom* var_chrom;
    // For adding substitutions:
    SubMutator subs;
    // For adding indels:
    IndelMutator indels;


    MutationSampler() : var_chrom(nullptr) {}

    MutationSampler(const MutationSampler& other)
        : var_chrom(other.var_chrom), subs(other.subs), indels(other.indels) {}

    MutationSampler& operator=(const MutationSampler& other) {
        var_chrom = other.var_chrom;
        subs = other.subs;
        indels = other.indels;
        return *this;
    }

    void new_chrom(VarChrom& var_chrom_, pcg64& eng) {
        var_chrom = &var_chrom_;
        subs.new_chrom(var_chrom_, eng);
        indels.new_chrom(var_chrom_);
        return;
    }
    // For writing to a file (used internally for testing):
    void new_chrom(VarChrom& var_chrom_, pcg64& eng, FileUncomp& file) {
        var_chrom = &var_chrom_;
        subs.new_chrom(var_chrom_, eng);
        subs.write_gammas(file);
        indels.new_chrom(var_chrom_);
        return;
    }

    /*
     Add mutations for a branch within a range.
     It also updates `end` for indels that occur in the range.
     (`end == begin` when chromosome region is of size zero bc `end` is non-inclusive)
     */
    void mutate(const double& b_len, pcg64& eng, const uint64& begin, uint64& end);

};












#endif
