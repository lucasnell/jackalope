#ifndef __GEMINO_MEVO_H
#define __GEMINO_MEVO_H


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution
#include <memory>  // smart pointers


#include "gemino_types.h" // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling



using namespace Rcpp;

namespace mevo {
    const std::string bases = "TCAG";
}


/*
 One sequence's Gamma values
 THIS IS TOO SIMPLE: INDELS CHANGE RATES OF SEQUENCES OUTSIDE THE INITIAL RANGE
 ... perhaps unless you use old positions rather than new ones to find Gamma values
 */
struct SeqGammas {
    std::vector<double> gamma_vals;  // values from a Gamma distribution
    uint k;                          // number of bp per value

    SeqGammas(const uint& k_, const RefSequence& rs,
              pcg32& eng, const double& alpha) : gamma_vals(), k(k_) {
        // Now fill gammas vector
        uint n_gammas = static_cast<uint>(std::ceil(
            static_cast<double>(rs.size()) / static_cast<double>(k_)));
        gamma_vals = std::vector<double>(n_gammas);
        std::gamma_distribution<double> distr(alpha, alpha);
        for (uint i = 0; i < n_gammas; i++) {
            gamma_vals[i] = distr(eng);
        }
    }

    double get_gamma(const uint& pos) {
        uint ind = pos / k;
        if (ind >= gamma_vals.size()) stop("ind too high for this SeqGammas object");
        return gamma_vals[ind];
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
 For table-sampling events depending on which nucleotide you start with.
 The `event_lengths` vector tells how long each event is.
 This field is 0 for substitions, < 0 for deletions, and > 0 for insertions.
 */
class MutationSampler {
public:

    std::vector<TableSampler> sampler;
    std::vector<sint> event_lengths;

    MutationSampler() : sampler(4), event_lengths() {
        base_inds = make_base_inds();
    }
    // copy constructor
    MutationSampler(const MutationSampler& other)
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

private:
    std::vector<uint8> base_inds;
};




class RateGetter {
public:
    RateGetter(const std::string& S_, const MutationRates& rates_)
    : S(S_), rates(rates_) {};
    inline double operator[](const uint& c) const {
        return rates[S[c]];
    }
    // Assignment operator
    RateGetter(const RateGetter& rhs) : S(rhs.S), rates(rhs.rates) {}

private:
    const std::string& S;
    const MutationRates& rates;
};






// MevoSampler combined objects for table-sampling event types and new nucleotides
class MevoSampler {
public:

    // For overall mutation rates by nucleotide:
    MutationRates rates;

    MevoSampler(const arma::mat& Q,
                const double& xi, const double& psi, const std::vector<double>& pis,
                arma::vec rel_insertion_rates, arma::vec rel_deletion_rates);
    MevoSampler(const MutationSampler& event_,
                const TableStringSampler<std::string>& nucleo_,
                const MutationRates& rates_)
        : rates(rates_), muts(event_), nts(nucleo_) {};

    /*
     Sample for mutation type based on nucleotide and rng engine
     */
    inline MutationInfo sample_muts(const char& c, pcg32& eng) const {
        return muts.sample(c, eng);
    }

    /*
     Create a new string of nucleotides (for insertions) of a given length and using
     an input rng engine
    */
    inline std::string new_nts(const uint& len, pcg32& eng) const {
        std::string str(len, 'x');
        nts.sample(str, eng);
        return str;
    }

private:
    // For sampling the type of mutation:
    MutationSampler muts;
    // For insertion sequences:
    TableStringSampler<std::string> nts;
};


#endif
