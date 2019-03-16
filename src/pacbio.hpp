#ifndef __GEMINO_PACBIO_H
#define __GEMINO_PACBIO_H




#include <RcppArmadillo.h>
#include <cmath>
#include <vector>  // vector class
#include <deque>  // deque class
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

// Used in PacBioReadLenSampler
#define SL_DEFAULT_S 0.200110276521
#define SL_DEFAULT_LOC -10075.4363813
#define SL_DEFAULT_SCALE 17922.611306
#define SL_DEFAULT_MIN_FRAG_LEN 50.0
// Used in PacBioPassSampler
#define SL_DEFAULT_CHI2_N1 0.00189237136
#define SL_DEFAULT_CHI2_N2 2.53944970
#define SL_DEFAULT_CHI2_N3 5500.0
#define SL_DEFAULT_CHI2_S1 0.01214
#define SL_DEFAULT_CHI2_S2 -5.12
#define SL_DEFAULT_CHI2_S3 675.0
#define SL_DEFAULT_CHI2_S4 48303.0732881
#define SL_DEFAULT_CHI2_S5 1.4691051212330266
#define SL_DEFAULT_MAX_PASSES 40
// Used in PacBioQualityError
#define SL_DEFAULT_SQRT_PARAMS1 0.5
#define SL_DEFAULT_SQRT_PARAMS2 0.2247
#define SL_DEFAULT_NORM_PARAMS1 0.0
#define SL_DEFAULT_NORM_PARAMS2 0.2
#define SL_DEFAULT_PROB_THRESH 0.2
#define SL_DEFAULT_PROB_INS 0.11
#define SL_DEFAULT_PROB_DEL 0.04
#define SL_DEFAULT_PROB_SUB 0.01





// Basic information to construct reads
struct PacBioReadConstrInfo {

    uint32 seq_ind;
    uint32 frag_len;
    uint32 frag_start;
    std::string read;
    std::string quals;
    uint32 read_seq_space;


    PacBioReadConstrInfo() {}
    PacBioReadConstrInfo(const PacBioReadConstrInfo& other)
        : seq_ind(other.seq_ind), frag_len(other.frag_len),
          frag_start(other.frag_start),
          read(other.read), quals(other.quals),
          read_seq_space(other.read_seq_space) {};
};








/*
 Sample for read lengths.

 If providing custom fragment lengths, make sure that the vector of their
 sampling weights is the same length as the vector of lengths.

 */
class PacBioReadLenSampler {

public:

    /* Initializers */
    // Default parameters from SimLoRD:
    PacBioReadLenSampler()
        : frag_lens(),
          sampler(),
          distr(std::log(SL_DEFAULT_SCALE), SL_DEFAULT_S),
          use_distr(true),
          min_frag_len(SL_DEFAULT_MIN_FRAG_LEN),
          loc(SL_DEFAULT_LOC) {};
    // Same but just modifying the min fragment length:
    PacBioReadLenSampler(const double& min_frag_len_)
        : frag_lens(),
          sampler(),
          distr(std::log(SL_DEFAULT_SCALE), SL_DEFAULT_S),
          use_distr(true),
          min_frag_len(std::ceil(min_frag_len_)),
          loc(SL_DEFAULT_LOC) {
        if (min_frag_len < 1) min_frag_len = 1;
    };
    // Using lognormal distribution with custom parameter values:
    PacBioReadLenSampler(const double& scale_,
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
    PacBioReadLenSampler(const std::vector<double>& probs,
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
    PacBioReadLenSampler(const PacBioReadLenSampler& other)
        : frag_lens(other.frag_lens),
          sampler(other.sampler),
          distr(other.distr),
          use_distr(other.use_distr),
          min_frag_len(other.min_frag_len),
          loc(other.loc) {};
    // Assignment operator
    PacBioReadLenSampler& operator=(const PacBioReadLenSampler& other) {
        frag_lens = other.frag_lens;
        sampler = other.sampler;
        distr = other.distr;
        use_distr = other.use_distr;
        min_frag_len = other.min_frag_len;
        loc = other.loc;
        return *this;
    }

    uint32 sample(pcg64& eng) {
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

private:

    std::vector<uint32> frag_lens;      // optional vector of possible fragment lengths
    TableSampler sampler;               // optional sampler that chooses from `frag_lens`
    std::lognormal_distribution<double> distr; // optional if using a distribution
    bool use_distr;                     // Whether to sample using `distr` field
    double min_frag_len;                // Minimum fragment length
    double loc;                         // Location parameter for lognormal distribution



};





/*
 Sample for number of passes.
 */
class PacBioPassSampler {

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


    void sample(uint32& split_pos,
                double& passes_left,
                double& passes_right,
                pcg64& eng,
                const double& read_length) {

        double passes, prop_left;

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

        passes = distr(eng);

        /*
         According to SimLoRD code, it's useful to not draw extreme outliers here
         because preventing those outliers "prevents outlier over the a/x boundary"
         */
        double outlier_threshold = R::qchisq(0.9925, n, 1, 0);
        while (passes > outlier_threshold) passes = distr(eng);
        // Now add scale and location parameters:
        passes *= s;
        passes += 1;
        // Make sure max pass isn't exceeded
        if (passes > max_passes) passes = max_passes;

        // Find point splitting between different # passes:
        double fraction, wholes;
        fraction = std::modf(passes, &wholes);


        if (static_cast<uint32>(wholes) % 2U == 0U) {
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

private:

    std::chi_squared_distribution<double> distr=std::chi_squared_distribution<double>(1);
    uint32 max_passes;
    std::vector<double> chi2_params_n;
    std::vector<double> chi2_params_s;

};





class PacBioQualityError {

public:

    // Initialize with default values:
    PacBioQualityError()
        : sqrt_params{SL_DEFAULT_SQRT_PARAMS1, SL_DEFAULT_SQRT_PARAMS2},
          norm_params{SL_DEFAULT_NORM_PARAMS1, SL_DEFAULT_NORM_PARAMS2},
          prob_thresh(SL_DEFAULT_PROB_THRESH),
          prob_ins(SL_DEFAULT_PROB_INS),
          prob_del(SL_DEFAULT_PROB_DEL),
          prob_subst(SL_DEFAULT_PROB_SUB),
          min_exp(calc_min_exp()) {};
    // Copy constructor
    PacBioQualityError(const PacBioQualityError& other)
        : sqrt_params(other.sqrt_params),
          norm_params(other.norm_params),
          prob_thresh(other.prob_thresh),
          prob_ins(other.prob_ins),
          prob_del(other.prob_del),
          prob_subst(other.prob_subst),
          min_exp(other.min_exp) {};
    // Assignment operator
    PacBioQualityError& operator=(const PacBioQualityError& other) {
        sqrt_params = other.sqrt_params;
        norm_params = other.norm_params;
        prob_thresh = other.prob_thresh;
        prob_ins = other.prob_ins;
        prob_del = other.prob_del;
        prob_subst = other.prob_subst;
        min_exp = other.min_exp;
        return *this;
    }
    // To limit the ways these can be changed:
    void change_sqrt_params(const double& sqrt_params1, const double& sqrt_params2) {
        sqrt_params = {sqrt_params1, sqrt_params2};
    }
    void change_norm_params(const double& norm_params1, const double& norm_params2) {
        norm_params = {norm_params1, norm_params2};
    }
    /*
     Change one or more error probabilities and the probability threshold.
     If you don't want one or more of them changed, just use a value outside
     the range [0,1).
     */
    void change_probs(const double& prob_thresh_,
                      const double& prob_ins_,
                      const double& prob_del_,
                      const double& prob_subst_) {
        if (prob_thresh_ >= 0 && prob_thresh_ < 1) prob_thresh = prob_thresh_;
        if (prob_ins_ >= 0 && prob_ins_ < 1) prob_ins = prob_ins_;
        if (prob_del_ >= 0 && prob_del_ < 1) prob_del = prob_del_;
        if (prob_subst_ >= 0 && prob_subst_ < 1) prob_subst = prob_subst_;
        min_exp = calc_min_exp();
    }


    void sample(pcg64& eng,
                char& qual_left,
                char& qual_right,
                std::deque<uint32>& insertions,
                std::deque<uint32>& deletions,
                std::deque<uint32>& substitutions,
                const uint32& read_length,
                const uint32& split_pos,
                const double& passes_left,
                const double& passes_right) {
        insertions.clear();
        deletions.clear();
        substitutions.clear();
        // Update sampling-error probabilities:
        modify_probs(eng, passes_left, passes_right);
        // Update qualities
        fill_quals(qual_left, qual_right);
        // Now iterate through and update insertions, deletions, and substitutions:
        uint32 current_length = 0;
        uint32 seq_pos = 0; // position on the sequence where events occur
        double u;
        std::vector<double>* cum_probs = &cum_probs_left;
        while (current_length < read_length) {
            if (current_length == split_pos) cum_probs = &cum_probs_right;
            u = runif_01(eng);
            if (u > cum_probs->at(2)) { // ------------ no errors
                current_length++;
            } else if (u < cum_probs->at(0)) { // ----- insertion
                // Don't add insertion if it would change read length
                if (current_length < (read_length - 1)) {
                    insertions.push_back(seq_pos);
                    current_length++;
                    if (current_length == split_pos) cum_probs = &cum_probs_right;
                }
                current_length++;
            } else if (u < cum_probs->at(1)) { // ----- deletion
                deletions.push_back(seq_pos);
            } else { // ------------------------------- substitution
                substitutions.push_back(seq_pos);
                current_length++;
            }
            seq_pos++;
        }
        return;
    }



private:

    // Values that change:
    std::vector<double> cum_probs_left = std::vector<double>(3);
    std::vector<double> cum_probs_right = std::vector<double>(3);
    // Values that do not change:
    std::vector<double> sqrt_params;
    std::vector<double> norm_params;
    double prob_thresh;
    double prob_ins;
    double prob_del;
    double prob_subst;
    double min_exp;   // always initialize last
    uint32 max_qual = 93;
    uint32 qual_start = static_cast<uint32>('!'); // quality of zero


    double calc_min_exp();

    inline double sigmoid(const double& x) {
        return 1 / (1 + std::pow(2, (-2.5 / 3 * x + 6.5 / 3)));
    }

    // Normal distribution truncated with lower threshold
    inline double trunc_norm(const double& lower_thresh,
                             pcg64& eng) {
        double rnd;
        double a_bar = (lower_thresh - norm_params[0]) / norm_params[1];

        /*
         The "near" method is preferred over the "far" method unless we're truncating
         all but the tail of the distribution:
         */
        if (lower_thresh < (norm_params[0] + 5 * norm_params[1])) {
            trunc_rnorm_near(rnd, a_bar, norm_params[0], norm_params[1], eng);
        } else {
            trunc_rnorm_far(rnd, a_bar, norm_params[0], norm_params[1], eng);
        }

        return rnd;
    }


    /*
    From SimLoRD:
    "Modify the given subread probabilities with an increase factor based on the
    number of passes. The increase is determined with a noisy sqare root function
    adapted with a sigmoidal factor. For the purpose of quality trimming the
    increase exponent is bounded with min_exp.
    Return the cumulative modfified probabilities for the left and right part of
    the read: (cum_probs_left, cum_probs_right) with (ins, ins+del, ins+del+subst)."
    */

    void modify_probs(pcg64& eng,
                      const double& passes_left,
                      const double& passes_right);


    void fill_quals(char& qual_left, char& qual_right) {
        uint32 tmp_l = std::round(-10.0 * std::log10(cum_probs_left.back()));
        uint32 tmp_r = std::round(-10.0 * std::log10(cum_probs_right.back()));
        if (tmp_l > max_qual) tmp_l = max_qual;
        if (tmp_r > max_qual) tmp_r = max_qual;
        qual_left = static_cast<char>(tmp_l + qual_start);
        qual_right = static_cast<char>(tmp_r + qual_start);
        return;
    }





};







/*
Template class to combine everything for Pac Bio sequencing of a single genome.
(We will need multiple of these objects to sequence a `VarSet` class.
See `PacBioVariants` class below.)

`T` should be `VarGenome` or `RefGenome`

*/
template <typename T>
class PacBioOneGenome {
public:

    /* __ Samplers __ */
    // Samples index for which genome-sequence to sequence
    AliasSampler seq_sampler;
    // Samples Pac Bio qualities and errors
    PacBioQualityError qe_sampler;
    // Samples numbers of passes over read:
    PacBioPassSampler pass_sampler;
    // Samples read lengths:
    PacBioReadLenSampler len_sampler;


    /* __ Info __ */
    std::vector<uint32> seq_lengths;    // genome-sequence lengths
    const T* sequences;                 // pointer to `const T`


    PacBioOneGenome() {};

    PacBioOneGenome(const PacBioOneGenome& other)
        : seq_sampler(other.seq_sampler),
          qe_sampler(other.qe_sampler),
          pass_sampler(other.pass_sampler),
          len_sampler(other.len_sampler),
          insertions(other.insertions),
          deletions(other.deletions),
          substitutions(other.substitutions),
          constr_info(other.constr_info) {};


    // Add one read string (with 4 lines: ID, sequence, "+", quality) to a FASTQ chunk
    void one_read(std::string& fastq_chunk,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info);

    /*
     Add information about a RefGenome or VarGenome object
     This is used when making multiple samplers that share most info except for
     that related to the sequence object.
     */
    void add_seq_info(const T& seq_object) {
        seq_lengths = seq_object.seq_sizes();
        sequences = &seq_object;

        construct_seqs();
    }


private:

    // To store indel locations, where each vector will be of length 2 if paired==true
    std::deque<uint32> insertions;
    std::deque<uint32> deletions;
    std::deque<uint32> substitutions;
    // Info to construct reads:
    PacBioReadConstrInfo constr_info;


    // Construct sequence-sampling probabilities:
    void construct_seqs() {
        std::vector<double> probs_;
        probs_.reserve(seq_lengths.size());
        for (uint i = 0; i < seq_lengths.size(); i++) {
            probs_.push_back(static_cast<double>(seq_lengths[i]));
        }
        seq_sampler = AliasSampler(probs_);
    }

};


typedef PacBioOneGenome<RefGenome> PacBioReference;
typedef PacBioOneGenome<VarGenome> PacBioOneVariant;











#endif
