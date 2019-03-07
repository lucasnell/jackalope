#ifndef __GEMINO_ALIAS_SAMPLER_H
#define __GEMINO_ALIAS_SAMPLER_H


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

#include "gemino_types.h" // integer types
#include "pcg.h"  // pcg seeding


using namespace Rcpp;





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



#endif
