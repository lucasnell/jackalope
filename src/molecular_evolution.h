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

    /*
     Below `q_tcag` is a length-4 vector of rates (q_i from Yang (2006)), for
     T, C, A, and G, respectively
     */
    MutationRates(const VarSequence& vs_, const std::vector<double>& q_tcag,
                  const SequenceGammas& gammas_)
        : vs(vs_), nt_rates(256, 0.0), gammas(gammas_) {
        for (uint i = 0; i < 4; i++) {
            uint j = mevo::bases[i];
            nt_rates[j] = q_tcag[i];
        }
    }
    MutationRates(const MutationRates& other)
        : vs(other.vs), nt_rates(other.nt_rates), gammas(other.gammas) {}

    // Using bracket operator to get the overall mutation rate at a location
    inline double operator[](const uint& pos) const {
        char c = vs.get_nt(pos);
        double r = nt_rates[c];
        r *= gammas[pos];
        return r;
    }

    /*
     The same as above, but for a range of positions.
     */
    inline double operator()(const uint& start, const uint& end) const {

        std::string seq;
        uint mut_;
        vs.set_seq_chunk(seq, start, end - start + 1, mut_);

        std::vector<double> gamma_vals = gammas(start, end);
        if (gamma_vals.size() != seq.size()) {
            stop("seq and gamma_vals sizes not matching in MutationRates::().");
        }

        double out = 0;
        for (uint i = 0; i < gamma_vals.size(); i++) {
            double r = nt_rates[seq[i]];
            out += r * gamma_vals[i];
        }

        return out;
    }

    // Return a single gamma value
    inline double gamma(const uint& pos) const {
        return gammas[pos];
    }


    // To get size of the variant sequence
    inline uint size() const noexcept {
        return vs.size();
    }

    // To return the overall rate for an entire sequence:
    double total_rate() const {

        double out = 0;

        if ((vs.size() - 1) != gammas.regions.back().end) {
            stop("gammas and vs sizes don't match inside MutationRates");
        }

        // If there are no mutations, this is pretty easy:
        if (vs.mutations.empty()) {

            uint end = vs.ref_seq.nucleos.size();
            if ((end - 1) != gammas.regions.back().end) {
                stop("gammas and vs ref sizes don't match inside MutationRates");
            }

            for (uint i = 0, idx = 0; idx < gammas.regions.size(); idx++) {
                double gamma = gammas.regions[idx].gamma;
                double tmp = 0;
                while (gammas.regions[idx].end >= i) {
                    tmp += nt_rates[vs.ref_seq.nucleos[i]];
                    i++;
                }
                out += (tmp * gamma);
            }

            return out;
        }

        /*
         std::vector<double> nt_rates;
         SequenceGammas gammas;
         const VarSequence& vs;
         */

        // Index to the first Mutation object
        uint mut_i = 0;
        // Index to the first gamma region
        uint gam_i = 0;
        double gamma = gammas.regions[gam_i].gamma;
        uint gamma_end = gammas.regions[gam_i].end;

        // Current position
        uint pos = 0;

        // Picking up any nucleotides before the first mutation
        for (; pos < vs.mutations[mut_i].new_pos; pos++) {
            if (pos > gamma_end) {
                gam_i++;
                gamma = gammas.regions[gam_i].gamma;
                gamma_end = gammas.regions[gam_i].end;
            }
            out += (nt_rates[vs.ref_seq[pos]] * gamma);
        }
        if (pos > gamma_end) {
            gam_i++;
            gamma = gammas.regions[gam_i].gamma;
            gamma_end = gammas.regions[gam_i].end;
        }


        // Now, for each subsequent mutation except the last, add all nucleotides
        // at or after its position but before the next one
        uint next_mut_i = mut_i + 1;
        while (next_mut_i < vs.mutations.size()) {
            while (pos < vs.mutations[next_mut_i].new_pos) {
                char c = vs.get_char_(pos, mut_i);
                out += nt_rates[c] * gamma;
                ++pos;
                if (pos > gamma_end) {
                    gam_i++;
                    gamma = gammas.regions[gam_i].gamma;
                    gamma_end = gammas.regions[gam_i].end;
                }
            }
            ++mut_i;
            ++next_mut_i;
        }

        // Now taking care of nucleotides after the last Mutation
        while (pos < vs.seq_size) {
            char c = vs.get_char_(pos, mut_i);
            out += nt_rates[c] * gamma;
            ++pos;
            if (pos > gamma_end) {
                gam_i++;
                gamma = gammas.regions[gam_i].gamma;
                gamma_end = gammas.regions[gam_i].end;
            }
        }

        return out;
    }

    // To update gamma boundaries when indels occur:
    inline void update_gamma_regions(const uint& pos, const sint& size_change) {
        gammas.update_gamma_regions(pos, size_change);
        return;
    }

};


/*
 This class uses the info above, plus a class and fxn from `weighted_reservoir.h` to
 do weighted reservoir sampling for a single location at which to put a mutation.
 The weights are based on the nucleotide and sequence region.

 Class `C` should be `ReservoirRates` or `ChunkReservoirRates`.

 If using `ReservoirRates`, this samples the entire sequence, which can be
 inefficient for large sequences.

 If using `ChunkReservoirRates`, this samples a set number of sequence locations
 at a time instead of the whole thing.
 It extracts a "chunk" for weighted sampling by using non-weighted sampling without
 replacement.
 For large sequences, this is much more efficient and, from my testing, produces
 similar results.
 */
template <template <typename> class C>
class OneSeqLocationSampler {

protected:

    C<MutationRates> rates;

public:

    OneSeqLocationSampler(const OneSeqLocationSampler<C>& other)
        : rates(other.rates) {}
    OneSeqLocationSampler(const MutationRates& mr, const uint chunk = 0)
        : rates(mr, chunk) {}

    inline uint sample(pcg32& eng) {
        return rates.sample(eng);
    }

    inline void update_gamma_regions(const uint& pos, const sint& size_change) {
        rates.update_gamma_regions(pos, size_change);
        return;
    }

};



/*
 Simplifying names and adding the following functionality:
 - return rate when using parentheses and a single uint argument
 - return a std::vector<double> of rates when using parentheses and two uint arguments
 - return a gamma value if using the `gamma()` method
 - return rate of the whole sequence with `total_rate()` method
 */
class LocationSampler: public OneSeqLocationSampler<ReservoirRates> {
public:

    // Constructor:
    LocationSampler(const MutationRates& mr) :
        OneSeqLocationSampler<ReservoirRates>(mr) {};

    /*
     Return a rate based on the actual indexing, rather than from `inds`.
     This was created to allow it to output a rate when you're adjusting the overall
     sequence rate during jumps across a phylogeny.
    */
    inline double operator()(const uint& pos) const {
        const MutationRates& mr(rates.res_rates);
        return mr[pos];
    }
    /*
    The same as above, but for a range of positions.
    */
    inline double operator()(const uint& start, const uint& end) const {
        const MutationRates& mr(rates.res_rates);
        return mr(start, end);
    }
    /*
     Return a single gamma value
     */
    inline double gamma(const uint& pos) const {
        const MutationRates& mr(rates.res_rates);
        return mr.gamma(pos);
    }
    /*
     Return the total rate for a VarSequence object
    */
    inline double total_rate() const {
        const MutationRates& mr(rates.res_rates);
        return mr.total_rate();
    }
};

class ChunkLocationSampler: public OneSeqLocationSampler<ChunkReservoirRates> {
public:

    // Constructor:
    ChunkLocationSampler(const MutationRates& mr, const uint chunk = 0)
        : OneSeqLocationSampler<ChunkReservoirRates>(mr, chunk) {}


    inline double operator()(const uint& pos) const {
        const MutationRates& mr(rates.res_rates.all_rates);
        return mr[pos];
    }

    inline double operator()(const uint& start, const uint& end) const {
        const MutationRates& mr(rates.res_rates.all_rates);
        return mr(start, end);
    }

    inline double gamma(const uint& pos) const {
        const MutationRates& mr(rates.res_rates.all_rates);
        return mr.gamma(pos);
    }

    inline double total_rate() const {
        const MutationRates& mr(rates.res_rates.all_rates);
        return mr.total_rate();
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
    // Initialize from an index and mut-lengths vector
    MutationInfo (const uint& ind, const std::vector<sint>& mut_lengths)
        : nucleo('\0'), length(0) {
        if (ind < 4) {
            nucleo = mevo::bases[ind];
        } else {
            length = mut_lengths[ind];
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
 The `mut_lengths` vector tells how long each mutation is.
 This field is 0 for substitions, < 0 for deletions, and > 0 for insertions.
 The `base_inds` field allows me to convert the characters 'T', 'C', 'A', or 'G'
 (cast to uints) into uints from 0 to 3.
 */
class MutationTypeSampler {

    std::vector<TableSampler> sampler;
    std::vector<sint> mut_lengths;
    std::vector<uint8> base_inds;

public:

    MutationTypeSampler() : sampler(4), mut_lengths(), base_inds() {
        base_inds = make_base_inds();
    }
    MutationTypeSampler(const std::vector<std::vector<double>>& probs,
                        const std::vector<sint>& mut_lengths_)
    : sampler(4), mut_lengths(mut_lengths_), base_inds() {
        base_inds = make_base_inds();
        if (probs.size() != 4) stop("probs must be size 4.");
        for (uint i = 0; i < 4; i++) {
            std::vector<double> probs_i = probs[i];
            /*
             Item `i` in `probs_i` needs to be manually converted to zero so
             it's not sampled.
             (Remember, we want the probability of each, *given that a mutation occurs*.
             Mutating into itself doesn't count.)
            */
            probs_i[i] = 0;
            // Now create and fill the `TableSampler`:
            sampler[i] = TableSampler(probs_i);
        }
    }
    // copy constructor
    MutationTypeSampler(const MutationTypeSampler& other)
        : sampler(other.sampler), mut_lengths(other.mut_lengths),
          base_inds(other.base_inds) {}
    // Assignment operator
    MutationTypeSampler& operator=(const MutationTypeSampler& other) {
        sampler = other.sampler;
        mut_lengths = other.mut_lengths;
        base_inds = other.base_inds;
        return *this;
    }

    /*
     Sample a mutation based on an input nucleotide.
     `c` gets cast to an uint, which is then input to `base_inds` to get the index
     from 0 to 3.
     */
    MutationInfo sample(const char& c, pcg32& eng) const {
        uint ind = sampler[base_inds[c]].sample(eng);
        MutationInfo mi(ind, mut_lengths);
        return mi;
    }

};









/*
 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================

 Combining samplers for location and for mutation type into a mutation sampler for
 a single sequence.

 =========================================================================================
 =========================================================================================
 =========================================================================================
 =========================================================================================
 */



/*
 OneSeqMutationSampler combines objects for sampling mutation types and new nucleotides for
 insertions.

 Class `C` should be `LocationSampler` or `ChunkLocationSampler`.
 */
template <class C>
class OneSeqMutationSampler {

    // VarSequence object reference to be manipulated
    VarSequence& vs;
    // For sampling the mutation location:
    C location;
    // For sampling the type of mutation:
    MutationTypeSampler type;
    // For new insertion sequences:
    TableStringSampler<std::string> insert;

    /*
     Sample for mutation location based on rates by sequence region and nucleotide.
     */
    inline uint sample_location(pcg32& eng) {
        return location.sample(eng);
    }

    /*
    Sample for mutation type based on nucleotide and rng engine
    */
    inline MutationInfo sample_type(const char& c, pcg32& eng) const {
        return type.sample(c, eng);
    }

    /*
    Create a new string of nucleotides (for insertions) of a given length and using
    an input rng engine
    */
    inline std::string new_nucleos(const uint& len, pcg32& eng) const {
        std::string str(len, 'x');
        insert.sample(str, eng);
        return str;
    }

public:

    // OneSeqMutationSampler(VarSequence& vs_) : vs(vs_) {}

    OneSeqMutationSampler(VarSequence& vs_,
                    const C& location_,
                    const MutationTypeSampler& type_,
                    const TableStringSampler<std::string>& insert_)
        : vs(vs_), location(location_),type(type_), insert(insert_) {}


    void mutate(pcg32& eng) {
        uint pos = sample_location(eng);
        char c = vs.get_nt(pos);
        MutationInfo m = sample_type(c, eng);
        if (m.length == 0) {
            vs.add_substitution(m.nucleo, pos);
        } else {
            if (m.length > 0) {
                std::string nts = new_nucleos(m.length, eng);
                vs.add_insertion(nts, pos);
            } else {
                if ((static_cast<sint>(pos) - m.length) > static_cast<sint>(vs.size())) {
                    m.length = static_cast<sint>(pos) - static_cast<sint>(vs.size());
                }
                uint del_size = std::abs(m.length);
                vs.add_deletion(del_size, pos);
            }
            // Update Gamma region bounds:
            location.update_gamma_regions(pos, m.length);
        }
        return;
    }

    // // Same as above, but it returns the change in the sequence rate that results
    // double mutate_rate_change(pcg32& eng) {
    //     uint pos = sample_location(eng);
    //     char c = vs.get_nt(pos);
    //     MutationInfo m = sample_type(c, eng);
    //     if (m.length == 0) {
    //         vs.add_substitution(m.nucleo, pos);
    //     } else {
    //         if (m.length > 0) {
    //             std::string nts = new_nucleos(m.length, eng);
    //             vs.add_insertion(nts, pos);
    //         } else {
    //             if ((static_cast<sint>(pos) - m.length) > static_cast<sint>(vs.size())) {
    //                 m.length = static_cast<sint>(pos) - static_cast<sint>(vs.size());
    //             }
    //             uint del_size = std::abs(m.length);
    //             vs.add_deletion(del_size, pos);
    //         }
    //         // Update Gamma region bounds:
    //         location.update_gamma_regions(pos, m.length);
    //     }
    //     return;
    // }
};


// Shortening these names
typedef OneSeqMutationSampler<LocationSampler> MutationSampler;
typedef OneSeqMutationSampler<ChunkLocationSampler> ChunkMutationSampler;




void fill_mut_prob_length_vectors(
        std::vector<std::vector<double>>& probs,
        std::vector<sint>& mut_lengths,
        const arma::mat& Q,
        const double& xi,
        const double& psi,
        const std::vector<double>& pi_tcag,
        arma::vec rel_insertion_rates,
        arma::vec rel_deletion_rates);


MutationSampler make_mutation_sampler(VarSequence& vs,
                                      const std::vector<std::vector<double>>& probs,
                                      const std::vector<sint>& mut_lengths,
                                      const std::vector<double>& pi_tcag,
                                      const arma::mat& gamma_mat);


ChunkMutationSampler make_mutation_sampler(VarSequence& vs,
                                           const std::vector<std::vector<double>>& probs,
                                           const std::vector<sint>& mut_lengths,
                                           const std::vector<double>& pi_tcag,
                                           const arma::mat& gamma_mat,
                                           const uint& chunk_size);



#endif
