#ifndef __JACKAL_ALIAS_SAMPLER_H
#define __JACKAL_ALIAS_SAMPLER_H


/*
 ********************************************************

 More numerically stable method of alias-sampling unsigned integers.
 From http://www.keithschwarz.com/darts-dice-coins/,
 section "Algorithm: Vose's Alias Method"

 ********************************************************
 */

#include <RcppArmadillo.h>
#include <vector>
#include <deque>
#include <string>
#include <pcg/pcg_random.hpp> // pcg prng

#include "jackal_types.h" // integer types
#include "pcg.h"  // pcg seeding
#include "util.h"  // str_stop


using namespace Rcpp;


namespace alias_sampler {

    const std::string bases = "TCAG";

}






class AliasSampler {
public:
    AliasSampler() : Prob(), Alias(), n(0) {};
    AliasSampler(const std::vector<double>& probs)
        : Prob(probs.size()), Alias(probs.size()), n(probs.size()) {

        arma::vec p(probs);
        p /= arma::accu(p);  // make sure they sum to 1
        p *= n;

        std::deque<uint32_t> Small;
        std::deque<uint32_t> Large;
        for (uint32_t i = 0; i < n; i++) {
            if (p(i) < 1) {
                Small.push_back(i);
            } else Large.push_back(i);
        }

        uint32_t l, g;
        while (!Small.empty() && !Large.empty()) {
            l = Small.front();
            Small.pop_front();
            g = Large.front();
            Large.pop_front();
            Prob[l] = p[l];
            Alias[l] = g;
            p(g) = (p(g) + p(l)) - 1;
            if (p(g) < 1) {
                Small.push_back(g);
            } else Large.push_back(g);
        }
        while (!Large.empty()) {
            g = Large.front();
            Large.pop_front();
            Prob[g] = 1;
        }
        while (!Small.empty()) {
            l = Small.front();
            Small.pop_front();
            Prob[l] = 1;
        }

    }
    // Copy constructor
    AliasSampler(const AliasSampler& other)
        : Prob(other.Prob), Alias(other.Alias), n(other.n) {}

    // Actual alias sampling
    inline uint32_t sample(pcg64& eng) const {
        // Fair dice roll from n-sided die
        uint32_t i = runif_01(eng) * n;
        // uniform in range (0,1)
        double u = runif_01(eng);
        if (u < Prob[i]) return(i);
        return Alias[i];
    };

private:
    std::vector<double> Prob;
    std::vector<uint32_t> Alias;
    uint32_t n;
};




/*
 Class template for table sampling a string, using an underlying AliasSampler object.
 `chars_in` should be the characters to sample from, `probs` the probabilities of
 sampling those characters.
 `T` can be `std::string` or `RefSequence`. Others may work, but are not guaranteed.
 */
template <typename T>
class AliasStringSampler {
public:

    T characters;

    AliasStringSampler(const T& chars_in, const std::vector<double>& probs)
        : characters(chars_in), uint_sampler(probs), n(probs.size()) {
        if (probs.size() != chars_in.size()) {
            str_stop({"For a AliasStringSampler construction, arguments probs and ",
                     "chars_in  must be same length."});
        }
    }
    AliasStringSampler() {}
    // copy constructor
    AliasStringSampler(const AliasStringSampler& other)
        : characters(other.characters), uint_sampler(other.uint_sampler),
          n(other.n) {}

    void sample(std::string& str, pcg64& eng) const {
        for (uint32 i = 0; i < str.size(); i++) {
            uint32 k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }
    char sample(pcg64& eng) const {
        uint32 k = uint_sampler.sample(eng);
        return characters[k];
    }

private:
    AliasSampler uint_sampler;
    uint32 n;
};





#endif
