#ifndef __JACKALOPE_MUTATOR_H
#define __JACKALOPE_MUTATOR_H


/*
 Combining samplers for location and for mutation type into a mutation sampler for
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
#include "mutator_location.h"  // mutator location sampling classes
#include "mutator_type.h"      // mutator type sampling classes
#include "util.h"  // str_stop



using namespace Rcpp;





/*
 MutationSampler combines objects for sampling mutation types and new
 nucleotides for insertions.
 */
class MutationSampler {

    /*
     Create a new string of nucleotides (for insertions) of a given length and using
     an input rng engine
     */
    inline std::string new_nucleos(const uint64& len, pcg64& eng) const {
        std::string str(len, 'x');
        insert.sample(str, eng);
        return str;
    }
    // Does most of the work for the mutate methods (all but location sampling)
    inline double mutate__(pcg64& eng, const uint64& start, sint64& end);

public:

    // VarChrom object pointer to be manipulated
    VarChrom* var_chrom;
    // For sampling the mutation location:
    LocationSampler location;
    // For sampling the type of mutation:
    MutationTypeSampler type;
    // For new insertion chromosomes:
    AliasStringSampler<std::string> insert;

    MutationSampler() : var_chrom(nullptr) {}

    MutationSampler(const MutationSampler& other)
        : var_chrom(other.var_chrom), location(other.location), type(other.type),
          insert(other.insert) {}

    MutationSampler& operator=(const MutationSampler& other) {
        if (other.var_chrom) var_chrom = other.var_chrom;
        location = other.location;
        type = other.type;
        insert = other.insert;
        return *this;
    }

    void new_chrom(VarChrom& vs_, const arma::mat& gamma_mat) {
        var_chrom = &vs_;
        location.new_chrom(vs_, gamma_mat);
        return;
    }

    // Add mutation and return the change in the chromosome rate that results
    double mutate(pcg64& eng);

    /*
     Overloaded for only mutating within a range.
     It also updates `end` if an indel occurs in the range.
     Make sure to keep checking for situation where `end < start` (i.e., chromosome section
     is empty).
     `// ***` mark difference between this and previous `mutate` versions
     */
    double mutate(pcg64& eng, const uint64& start, sint64& end);

};












#endif
