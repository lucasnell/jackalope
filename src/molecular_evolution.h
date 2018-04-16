#ifndef __GEMINO_MEVO_H
#define __GEMINO_MEVO_H


#include <RcppArmadillo.h>
#include <sitmo.h>    // sitmo
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // gamma_distribution
#include <algorithm>  // lower_bound, random_shuffle
#include <cmath>  // std::exp, std::log
#include <numeric>  // accumulate
#include <unordered_map>  // unordered_map
#include <queue>  // priority_queue
#ifdef _OPENMP
#include <omp.h>  // omp
#endif



#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "alias.h"    // Alias sampling


using namespace Rcpp;

namespace mevo {
    const std::string bases = "ACGT";
    std::unordered_map<char, uint> base_inds = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    const double sitmo_max = static_cast<double>(sitmo::prng_engine::max()) + 1.0;
}


/*
 One sequence's Gamma values
 THIS IS TOO SIMPLE: INDELS CHANGE RATES OF SEQUENCES OUTSIDE THE INITIAL RANGE
 */
struct SeqGammas {
    std::vector<double> gamma_vals;  // values from a Gamma distribution
    uint k;                          // number of bp per value

    SeqGammas(const uint& k_, const RefSequence& rs,
              sitmo::prng_engine& eng, const double& alpha) : gamma_vals(), k(k_) {
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




struct QMaps {

    std::unordered_map<char, std::vector<double>> Q;
    std::unordered_map<char, std::vector<double>> M;
    std::unordered_map<char, double> q;
    std::unordered_map<char, double> w;  // same as above, but sums to 1

    QMaps() {
        Q = {{'A', std::vector<double>(0)}, {'C', std::vector<double>(0)},
             {'G', std::vector<double>(0)}, {'T', std::vector<double>(0)}};
        M = {{'A', std::vector<double>(0)}, {'C', std::vector<double>(0)},
             {'G', std::vector<double>(0)}, {'T', std::vector<double>(0)}};
        q = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
        w = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    };
    QMaps(uint N) {
        Q = {{'A', std::vector<double>(N, 0.0)}, {'C', std::vector<double>(N, 0.0)},
             {'G', std::vector<double>(N, 0.0)}, {'T', std::vector<double>(N, 0.0)}};
        M = {{'A', std::vector<double>(N, 0.0)}, {'C', std::vector<double>(N, 0.0)},
             {'G', std::vector<double>(N, 0.0)}, {'T', std::vector<double>(N, 0.0)}};
        q = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
        w = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    };

};




// struct to store info on a nucleotide's key value and position for weighted sampling
struct NucleoKeyPos {

    double key;
    uint pos;

    NucleoKeyPos() : key(0.0), pos(0) {};
    NucleoKeyPos(const double& _key, const uint& _pos) : key(_key), pos(_pos) {};

};


/*
 For alias-sampling events depending on which nucleotide you start with
 The event_types vector indicates which type of event a location in one of
     the unordered_map vectors is: 0:3 for substitution to A, C, G, and T,
     respectively, 4 for insertion, and 5 for deletion.
     It should be the same length as each vector for a given nucleotide in the
     `F` or `L` fields.
 The event_lengths vector should be the coincide with items in the event_types
     vector, but it tells how long each event is. For substitutions, this field
     is ignored, but a filler value should be provided so the event_lengths
     and event_types vectors are the same length.
 */
struct EventAliasTables {

    std::unordered_map<char, std::vector<double>> F;
    std::unordered_map<char, std::vector<uint>> L;
    std::vector<uint> event_types;
    std::vector<uint> event_lengths;

    EventAliasTables() {
        F = {{'A', std::vector<double>(0)}, {'C', std::vector<double>(0)},
             {'G', std::vector<double>(0)}, {'T', std::vector<double>(0)}};
        L = {{'A', std::vector<uint>(0)}, {'C', std::vector<uint>(0)},
             {'G', std::vector<uint>(0)}, {'T', std::vector<uint>(0)}};
        event_types = std::vector<uint>(0);
        event_lengths = std::vector<uint>(0);
    };
    EventAliasTables(uint N) {
        F = {{'A', std::vector<double>(N, 0.0)}, {'C', std::vector<double>(N, 0.0)},
             {'G', std::vector<double>(N, 0.0)}, {'T', std::vector<double>(N, 0.0)}};
        L = {{'A', std::vector<uint>(N, 0)}, {'C', std::vector<uint>(N, 0)},
             {'G', std::vector<uint>(N, 0)}, {'T', std::vector<uint>(N, 0)}};
        event_types = std::vector<uint>(N, 0);
        event_lengths = std::vector<uint>(N, 0);
    }
};

// For alias-sampling nucleotides for new sequences (e.g. insertions)
struct NucleoAliasTables {

    std::vector<double> F;
    std::vector<uint> L;

    NucleoAliasTables() {
        F = std::vector<double>(4);
        L = std::vector<uint>(4);
    };
};

// Alias tables for both event types and nucleotide types
struct Alias {
    EventAliasTables event;
    NucleoAliasTables nucleo;

    Alias(EventAliasTables event_, NucleoAliasTables nucleo_)
        : event(event_), nucleo(nucleo_) {};
    Alias(uint N) : event(N), nucleo() {};
    Alias() : event(), nucleo() {};

};


#endif
