#ifdef __GEMINO_ABSOLUTELY_NEVER_COMPILE_H

/*
 Above line, plus extra `#endif` below should keep the rest of this file
 from being compiled.
 */



#ifndef __GEMINO_LONGREADS_H
#define __GEMINO_LONGREADS_H




#include <RcppArmadillo.h>
#include <cmath>
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng
#include <string>  // string class
#include <random>  // distributions

#include "gemino_types.hpp"  // uint32
#include "seq_classes_ref.hpp"  // Ref* classes
#include "seq_classes_var.hpp"  // Var* classes
#include "pcg.hpp"  // runif_01
#include "table_sampler.hpp"  // TableSampler
#include "alias_sampler.hpp"  // AliasSampler
#include "util.hpp"  // clear_memory
#include "str_manip.hpp"  // rev_comp
#include "sequencer.hpp"  // SequenceIdentifierInfo

using namespace Rcpp;



// Defaults from SimLoRD:

#define SL_DEFAULT_S = 0.200110276521
#define SL_DEFAULT_LOC = -10075.4363813
#define SL_DEFAULT_SCALE = 17922.611306
#define SL_DEFAULT_MIN_FRAG_LEN = 50
#define SL_DEFAULT_CHI2_N1 0.00189237136
#define SL_DEFAULT_CHI2_N2 2.53944970
#define SL_DEFAULT_CHI2_N3 5500
#define SL_DEFAULT_CHI2_S1 0.01214
#define SL_DEFAULT_CHI2_S2 -5.12
#define SL_DEFAULT_CHI2_S3 675
#define SL_DEFAULT_CHI2_S4 48303.0732881
#define SL_DEFAULT_CHI2_S5 1.4691051212330266
#define SL_DEFAULT_MAX_PASSES 40

// Not yet implemented:
#define SL_DEFAULT_SQRT_PARAMS1 0.5
#define SL_DEFAULT_SQRT_PARAMS2 0.2247
#define SL_DEFAULT_NORM_PARAMS1 0
#define SL_DEFAULT_NORM_PARAMS2 0.2
#define SL_DEFAULT_PROB_THRESH 0.2
#define SL_DEFAULT_PROB_INS 0.11
#define SL_DEFAULT_PROB_DEL 0.04
#define SL_DEFAULT_PROB_SUB 0.01






/*
 Sample for fragment lengths.

 If providing custom fragment lengths, make sure that the vector of their
 sampling weights is the same length as the vector of lengths.

 */
class PacBioFragLenSampler {


    std::vector<uint32> frag_lens;      // optional vector of possible fragment lengths
    TableSampler sampler;               // optional sampler that chooses from `frag_lens`
    std::lognormal_distribution<double> distr; // optional if using a distribution
    bool use_distr;                     // Whether to sample using `distr` field
    double min_frag_len;                // Minimum fragment length
    double loc;                         // Location parameter for lognormal distribution

public:

    /* Initializers */
    // Default parameters from SimLoRD:
    PacBioFragLenSampler()
        : frag_lens(),
          sampler(),
          distr(std::log(SL_DEFAULT_SCALE), SL_DEFAULT_S),
          use_distr(true),
          min_frag_len(SL_DEFAULT_MIN_FRAG_LEN),
          loc(SL_DEFAULT_LOC) {};
    // Same but just modifying the min fragment length:
    PacBioFragLenSampler(const double& min_frag_len_)
        : frag_lens(),
          sampler(),
          distr(std::log(SL_DEFAULT_SCALE), SL_DEFAULT_S),
          use_distr(true),
          min_frag_len(std::ceil(min_frag_len_)),
          loc(SL_DEFAULT_LOC) {
        if (min_frag_len < 1) min_frag_len = 1;
    };
    // Using lognormal distribution with custom parameter values:
    PacBioFragLenSampler(const double& scale_,
                         const double& sigma_,
                         const double& loc_,
                         const double& min_frag_len_)
        : frag_lens(),
          sampler(),
          distr(std::log(scale_), sigma_),
          use_distr(true),
          min_frag_len(std::ceil(min_frag_len_)),
          loc(loc_) {
        if (min_frag_len < 1) min_frag_len = 1;
    };
    // Using a vector of fragment lengths, each with a sampling probability:
    PacBioFragLenSampler(const std::vector<double>& probs,
                         const std::vector<uint32>& frag_lens_)
        : frag_lens(frag_lens_),
          sampler(probs),
          distr(),
          use_distr(false),
          min_frag_len(),
          loc() {
        if (probs.size() != frag_lens_.size()) {
            stop("Probability and fragment lengths vector should be the same length.");
        }
    };
    // Copy constructor
    PacBioFragLenSampler(const PacBioFragLenSampler& other)
        : frag_lens(other.frag_lens),
          sampler(other.sampler),
          distr(other.distr),
          use_distr(other.use_distr),
          min_frag_len(other.min_frag_len),
          loc(other.loc) {};
    // Assignment operator
    PacBioFragLenSampler& operator=(const PacBioFragLenSampler& other) {
        frag_lens = other.frag_lens;
        sampler = other.sampler;
        distr = other.distr;
        use_distr = other.use_distr;
        min_frag_len = other.min_frag_len;
        loc = other.loc;
        return *this;
    }

    uint32 sample(pcg64& eng) const {
        uint32 len_;
        if (use_distr) {
            double rnd = distr(eng) + loc;
            uint32 iters = 0; // to make sure it doesn't run many times
            // Rejection sampling to keep it above minimum:
            while (rnd < min_frag_len && iters < 10) {
                rnd = distr(eng) + loc;
                iters++;
            }
            // Give up if it keeps returning craziness:
            if (rnd < min_frag_len) rnd = min_frag_len;
            len_ = static_cast<uint32>(rnd);
        } else {
            uint64 ind = sampler.sample(eng);
            uint32 len_ = frag_lens[ind];
        }
        return len_;
    }

};





/*
 Sample for number of passes.
 */
class PacBioPassSampler {

    std::chi_squared_distribution<double> distr(1);
    uint32 max_passes;
    std::vector<double> chi2_params_n;
    std::vector<double> chi2_params_s;

public:

    /* Initializers */
    // All defaults:
    PacBioPassSampler()
        : max_passes(SL_DEFAULT_MAX_PASSES),
          chi2_params_n{SL_DEFAULT_CHI2_N1, SL_DEFAULT_CHI2_N2, SL_DEFAULT_CHI2_N3},
          chi2_params_s{SL_DEFAULT_CHI2_S1, SL_DEFAULT_CHI2_S2, SL_DEFAULT_CHI2_S3,
                        SL_DEFAULT_CHI2_S4, SL_DEFAULT_CHI2_S5} {};
    // Specifying for max passes only:
    PacBioPassSampler(const uint32& max_passes_)
        : max_passes(max_passes_),
          chi2_params_n{SL_DEFAULT_CHI2_N1, SL_DEFAULT_CHI2_N2, SL_DEFAULT_CHI2_N3},
          chi2_params_s{SL_DEFAULT_CHI2_S1, SL_DEFAULT_CHI2_S2, SL_DEFAULT_CHI2_S3,
                        SL_DEFAULT_CHI2_S4, SL_DEFAULT_CHI2_S5} {};
    // Copy constructor
    PacBioPassSampler(const PacBioPassSampler& other)
        : max_passes(other.max_passes),
          chi2_params_n(other.chi2_params_n),
          chi2_params_s(other.chi2_params_s) {};
    // Assignment operator
    PacBioPassSampler& operator=(const PacBioPassSampler& other) {
        max_passes = other.max_passes;
        chi2_params_n = other.chi2_params_n;
        chi2_params_s = other.chi2_params_s;
        return *this;
    }

    // Function to set just parameters for `n`
    void set_n(const double& chi2_n1,
               const double& chi2_n2,
               const double& chi2_n3) {
        chi2_params_n = {chi2_n1, chi2_n2, chi2_n3};
        return;
    }

    // Function to set just parameters for `s`
    void set_s(const double& chi2_s1,
               const double& chi2_s2,
               const double& chi2_s3,
               const double& chi2_s4,
               const double& chi2_s5) {
        chi2_params_s = {chi2_s1, chi2_s2, chi2_s3, chi2_s4, chi2_s5};
        return;
    }


    void sample(double& passes,
                uint32& split_pos,
                uint32& passes_left,
                uint32& passes_right,
                double& prop_left,
                pcg64& eng,
                const double& read_length) {

        double n = chi2_params_n[0] * std::min(read_length, chi2_params_n[2]) +
            chi2_params_n[1];
        if (n < 0.001) n = 0.001;

        double s;
        if (read_length <= chi2_params_s[2]) {
            s = chi2_params_s[0] * read_length - chi2_params_s[1];
            if (s < 0.001) s = 0.001;
        } else {
            s = chi2_params_s[3] / std::pow(read_length, chi2_params_s[4]);
        }

        // Reset distr for new `n`:
        distr.param(std::chi_squared_distribution<double>::param_type(n));

        double passes = distr(eng);

        /*
         According to SimLoRD code, it's useful to not draw extreme outliers here
         because preventing those outliers "prevents outlier over the a/x boundary"
         */
        double outlier_threshold = R::qchisq(0.9925, n, 1, 0);
        while (passes > outlier_threshold) pass = distr(eng);
        // Now add scale and location parameters:
        passes *= s;
        passes += 1;
        // Make sure max pass isn't exceeded
        if (passes > max_passes) passes = max_passes;

        // Find point splitting between different # passes:
        double fraction, wholes;
        fraction = std::modf(passes, &wholes);


        if (wholes % 2 == 0) {
            prop_left = fraction;
            split_pos = std::round(static_cast<double>(read_length) * prop_left);
            passes_left = std::ceil(passes);
            passes_right = std::floor(passes);
        } else {
            prop_left = 1 - fraction;
            split_pos = std::round(static_cast<double>(read_length) * prop_left);
            passes_left = std::floor(passes);
            passes_right = std::ceil(passes);
        }

        return;

    }

};




#endif

#endif
