#ifndef __GEMINO_MEVO_H
#define __GEMINO_MEVO_H


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution


#include "gemino_types.h" // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "mevo_gammas.h"  // SequenceGammas class
#include "weighted_reservoir.h"  // weighted_reservoir_* functions



using namespace Rcpp;

namespace mevo {
    const std::string bases = "TCAG";
}


/*
 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================

 Choosing mutation locations based on overall mutation rates that vary by
 (i) nucleotide and (ii) sequence region

 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================
 */






/*
 Stores info on the overall mutation rates (including Gammas) for each nucleotide.

 This class is used to do weighted reservoir sampling for where a mutation should occur.
 Ultimately this class allows you to use a bracket operator to get a rate for a
 given location in the variant genome.

 For nucleotide rates (not including Gammas), Ns are set to 0 bc we don't want to
 process these.
 Input char objects are cast to uint which provide the indices.
 T, C, A, G, and N should never be higher than 84, so will be safe.
 In case someone inputs other characters accidentally, I've set the length to 256,
 which should work for all 8-bit character values.
 The memory overhead should be pretty minimal.
 */

class MutationRates {

private:
    const VarSequence& vs;
    std::vector<double> nt_rates;
    SequenceGammas gammas;

public:
    MutationRates(const VarSequence& vs_, const std::vector<double>& rates_,
                  const uint& gamma_size_,
                  pcg32& eng, const double& alpha)
        : vs(vs_), nt_rates(256, 0.0), gammas(vs_, gamma_size_, eng, alpha) {
        for (uint i = 0; i < 4; i++) {
            uint j = mevo::bases[i];
            nt_rates[j] = rates_[i];
        }
    }
    MutationRates(const VarSequence& vs_, const std::vector<double>& rates_,
                  arma::mat gamma_mat)
        : vs(vs_), nt_rates(256, 0.0), gammas(gamma_mat) {
        for (uint i = 0; i < 4; i++) {
            uint j = mevo::bases[i];
            nt_rates[j] = rates_[i];
        }
    }
    // Assignment operator
    MutationRates(const MutationRates& rhs)
        : vs(rhs.vs), nt_rates(rhs.nt_rates), gammas(rhs.gammas) {}

    // Using bracket operator to get the overall mutation rate at a location
    inline double operator[](const uint& new_pos) const {
        char c = vs.get_nt(new_pos);
        double r = nt_rates[c];
        r *= gammas[new_pos];
        return r;
    }
    // To get size of the variant sequence
    inline uint size() const noexcept {
        return vs.size();
    }

};


/*
 This class uses the info above, plus a class and fxn from `weighted_reservoir.h` to
 do weighted reservoir sampling for a single location at which to put a mutation.
 This sampling is done using the entire sequence, which can be inefficient for large
 sequences.
 The weights are based on the nucleotide and sequence region.
 */
class LocationSampler {

    ReservoirRates<MutationRates> rates;

    LocationSampler(const VarSequence& vs_, const std::vector<double>& rates_,
                    arma::mat gamma_mat)
        : rates(MutationRates(vs_, rates_, gamma_mat)) {}
    LocationSampler(const VarSequence& vs_, const std::vector<double>& rates_,
                    const uint& gamma_size_,
                    pcg32& eng, const double& alpha)
        : rates(MutationRates(vs_, rates_, gamma_size_, eng, alpha)) {}

    inline uint sample(pcg32& eng) {
        return weighted_reservoir_<MutationRates>(rates, eng);
    }

};


/*
 This is the same as above, but it only does weighted sampling on a set number of
 sequence locations at a time.
 It extracts a "chunk" for weighted sampling by using non-weighted sampling without
 replacement.
 For large sequences, this is much more efficient and, from my testing, produces
 similar results.
 */
class ChunkLocationSampler {

    ChunkReservoirRates<MutationRates> rates;

    ChunkLocationSampler(const uint& chunk_size,
                         const VarSequence& vs_, const std::vector<double>& rates_,
                         arma::mat gamma_mat)
        : rates(MutationRates(vs_, rates_, gamma_mat), chunk_size) {}

    ChunkLocationSampler(const uint& chunk_size,
                         const VarSequence& vs_, const std::vector<double>& rates_,
                         const uint& gamma_size_,
                         pcg32& eng, const double& alpha)
        : rates(MutationRates(vs_, rates_, gamma_size_, eng, alpha), chunk_size) {}

    inline uint sample(pcg32& eng) {
        return weighted_reservoir_chunk_<MutationRates>(rates, eng);
    }

};







/*
 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================

 Choosing mutation type based on the starting nucleotide

 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================
 */



/*
 For a single mutation's info
 */
struct MutationInfo {
    char nucleo;
    sint length;
    // Initialize from an index and event-lengths vector
    MutationInfo (const uint& ind, const std::vector<sint>& event_lengths)
        : nucleo('\0'), length(0) {
        if (ind < 4) {
            nucleo = mevo::bases[ind];
        } else {
            length = event_lengths[ind];
        }
    }
};



/*
 For constructors, this creates a vector of indices for each char in "TCAG" (0 to 3).
 Because char objects can be easily cast to uints, I can input a char from a sequence
 and get out an index to which TableSampler object to sample from.
 This way is much faster than using an unordered_map.
 Using 8-bit uints bc the char should never be >= 256.
 It's only of length 85 (versus 256 in MutationRates class) because characters other
 than T, C, A, or G have rates hard-coded to zero and should never be chosen.
 */
inline std::vector<uint8> make_base_inds() {
    std::vector<uint8> base_inds(85);
    uint8 i = 0;
    for (const char& c : mevo::bases) {
        base_inds[c] = i;
        i++;
    }
    return base_inds;
}

/*
 For table-sampling mutation types depending on which nucleotide you start with.
 The `event_lengths` vector tells how long each event is.
 This field is 0 for substitions, < 0 for deletions, and > 0 for insertions.
 The `base_inds` field allows me to convert the characters 'T', 'C', 'A', 'G', or 'N'
 (cast to uints) into uints from 0 to 3.
 */
class MutationTypeSampler {

    std::vector<uint8> base_inds;
    std::vector<TableSampler> sampler;
    std::vector<sint> event_lengths;

public:

    MutationTypeSampler() : sampler(4), event_lengths() {
        base_inds = make_base_inds();
    }
    // copy constructor
    MutationTypeSampler(const MutationTypeSampler& other)
        : sampler(other.sampler), event_lengths(other.event_lengths),
          base_inds(other.base_inds) {}

    /*
     Sample an event based on an input nucleotide.
     `c` gets cast to an uint, which is then input to `base_inds` to get the index
     from 0 to 3.
     */
    MutationInfo sample(const char& c, pcg32& eng) const {
        uint ind = sampler[base_inds[c]].sample(eng);
        MutationInfo mi(ind, event_lengths);
        return mi;
    }

};






/*
 MutationSampler combines objects for sampling event types and new nucleotides for
 insertions.
 */
class MutationSampler {

private:

    // For sampling the type of mutation:
    MutationTypeSampler types;
    // For insertion sequences:
    TableStringSampler<std::string> nucleos;

public:

    // For overall mutation rates by nucleotide:
    MutationRates rates;

    MutationSampler(const arma::mat& Q,
                const double& xi, const double& psi, const std::vector<double>& pis,
                arma::vec rel_insertion_rates, arma::vec rel_deletion_rates);
    MutationSampler(const MutationTypeSampler& event_,
                const TableStringSampler<std::string>& nucleo_,
                const MutationRates& rates_)
        : rates(rates_), types(event_), nucleos(nucleo_) {};

    /*
     Sample for mutation type based on nucleotide and rng engine
     */
    inline MutationInfo sample_types(const char& c, pcg32& eng) const {
        return types.sample(c, eng);
    }

    /*
     Create a new string of nucleotides (for insertions) of a given length and using
     an input rng engine
    */
    inline std::string new_nucleos(const uint& len, pcg32& eng) const {
        std::string str(len, 'x');
        nucleos.sample(str, eng);
        return str;
    }

};


#endif
