#ifndef __JACKALOPE_ALIAS_SAMPLER_H
#define __JACKALOPE_ALIAS_SAMPLER_H


#include "jackalope_config.h" // controls debugging and diagnostics output

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

#include "jackalope_types.h" // integer types
#include "pcg.h"  // pcg seeding
#include "util.h"  // str_stop


using namespace Rcpp;







class AliasSampler {
public:
    AliasSampler() : Prob(), Alias(), n(0) {};
    AliasSampler(const std::vector<double>& probs)
        : Prob(probs.size()), Alias(probs.size()), n(probs.size()) {
        arma::rowvec p(probs);
        construct(p);
    }
    AliasSampler(arma::rowvec probs)
        : Prob(probs.n_elem), Alias(probs.n_elem), n(probs.n_elem) {
        construct(probs);
    }
    // Copy constructor
    AliasSampler(const AliasSampler& other)
        : Prob(other.Prob), Alias(other.Alias), n(other.n) {}

    // Actual alias sampling
    inline uint64 sample(pcg64& eng) const {
        // Fair dice roll from n-sided die
        uint64 i = runif_01(eng) * n;
        // uniform in range (0,1)
        double u = runif_01(eng);
        if (u < Prob[i]) return(i);
        return Alias[i];
    };

private:
    std::vector<double> Prob;
    std::vector<uint64> Alias;
    uint64 n;


    void construct(arma::rowvec& p) {

        p /= arma::accu(p);  // make sure they sum to 1
        p *= n;

        std::deque<uint64> Small;
        std::deque<uint64> Large;
        for (uint64 i = 0; i < n; i++) {
            if (p(i) < 1) {
                Small.push_back(i);
            } else Large.push_back(i);
        }

        uint64 l, g;
        while (!Small.empty() && !Large.empty()) {
            l = Small.front();
            Small.pop_front();
            g = Large.front();
            Large.pop_front();
            Prob[l] = p(l);
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

        return;
    }
};




/*
 Class template for table sampling a string, using an underlying AliasSampler object.
 `chars_in` should be the characters to sample from, `probs` the probabilities of
 sampling those characters.
 `T` can be `std::string` or `RefChrom`. Others may work, but are not guaranteed.
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
        for (uint64 i = 0; i < str.size(); i++) {
            uint64 k = uint_sampler.sample(eng);
            str[i] = characters[k];
        }
        return;
    }
    char sample(pcg64& eng) const {
        uint64 k = uint_sampler.sample(eng);
        return characters[k];
    }

private:
    AliasSampler uint_sampler;
    uint64 n;
};





#endif
