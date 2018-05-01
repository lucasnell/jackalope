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
#include "gammas.h"  // SequenceGammas class



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
 Chooses a chunk while controlling for the last chunk's size being < `chunk_size`.
 */
class ChunkChooser {

private:

    uint n_chunks;    // Number of chunks to choose from
    uint chunk_size;  // Size of each chunk
    double p_retain;  // Pr(retaining last chunk | last chunk was sampled initially):

public:

    ChunkChooser(const uint& chunk_size_, const uint& overall_size)
        : n_chunks(), chunk_size(chunk_size_), p_retain(overall_size) {

        n_chunks = std::ceil(static_cast<double>(overall_size) /
            static_cast<double>(chunk_size));

        p_retain -= static_cast<double>((n_chunks - 1) * chunk_size);
        p_retain /= overall_size;
        p_retain *= static_cast<double>(n_chunks);

    }

    uint choose(pcg32& eng) {

        if (n_chunks == 1) return 0;

        double u = runif_01(eng);
        uint chunk = u * n_chunks;
        if (chunk == n_chunks - 1) {
            u = runif_01(eng);
            if (u > p_retain) {
                u = runif_01(eng);
                chunk = u * (n_chunks - 1);
            }
        }
        return chunk;
    }

    /*
     Reset sizes after you add an insertion or deletion
     */
    void reset_size(const uint& new_size) {

        n_chunks = std::ceil(static_cast<double>(new_size) /
            static_cast<double>(chunk_size));

        p_retain = static_cast<double>(new_size);
        p_retain -= static_cast<double>((n_chunks - 1) * chunk_size);
        p_retain /= static_cast<double>(new_size);
        p_retain *= static_cast<double>(n_chunks);
    }
};






/*
 Stores info on the overall mutation rates for each nucleotide.
 Ns are set to 0 bc we don't want to process these.
 Input char objects are cast to uint which provide the indices.
 T, C, A, G, and N should never be higher than 84, so will be safe.
 If you're worried about other characters accidentally being input to it, you can
 set `rates` to size 256.
 */
class MutationRates {

public:

    std::vector<double> rates;

    MutationRates() : rates(85, 0.0) {}

    MutationRates(const std::vector<double>& rates_) : rates(85, 0.0) {
        for (uint i = 0; i < 4; i++) {
            uint j = mevo::bases[i];
            rates[j] = rates_[i];
        }
    }

    /*
     Return rates when square brackets are used
     `c` is intended to be a char converted to an integer.
     */
    double operator[](const uint& i) const {
        return rates[i];
    }
    double& operator[](const uint& i) {
        return rates[i];
    }
};




/*
 This class allows me to use a template in `molecular_evolution.cpp` to
 do weighted reservoir sampling.
 RateGetter combines references to a string and to a MutationRates object, and
 ultimately allows you to use a bracket operator to get a rate for a given location
 in a sequence.
 */
class RateGetter {

private:
    const std::string& S;
    const MutationRates& rates;

public:
    RateGetter(const std::string& S_, const MutationRates& rates_)
        : S(S_), rates(rates_) {};
    inline double operator[](const uint& idx) const {
        return rates[S[idx]];
    }
    // Assignment operator
    RateGetter(const RateGetter& rhs) : S(rhs.S), rates(rhs.rates) {}

};



// /*
//  This class allows me to use a template in `molecular_evolution.cpp` to
//  do weighted reservoir sampling.
//  RateGetter combines references to a string and to a MutationRates object, and
//  ultimately allows you to use a bracket operator to get a rate for a given location
//  in a sequence.
//  */
// class GammaGetter {
//
// private:
//     const MutationRates& rates;
//     const std::vector<double>& gammas;
//
// public:
//     RateGetter(const std::string& S_, const MutationRates& rates_)
//         : S(S_), rates(rates_) {};
//     inline double operator[](const uint& idx) const {
//         return rates[S[idx]];
//     }
//     // Assignment operator
//     RateGetter(const RateGetter& rhs) : S(rhs.S), rates(rhs.rates) {}
//
// };




// VarSequence








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

private:

    std::vector<uint8> base_inds;

public:

    std::vector<TableSampler> sampler;
    std::vector<sint> event_lengths;

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
