#ifndef __JACKALOPE_PACBIO_H
#define __JACKALOPE_PACBIO_H



#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>  // vector class
#include <deque>  // deque class
#include <pcg/pcg_random.hpp> // pcg prng
#include <string>  // string class
#include <random>  // distributions
#include <fstream> // for writing FASTQ files
#include "zlib.h"  // for writing to compressed FASTQ
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // uint64
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // runif_01
#include "alias_sampler.h"  // AliasSampler
#include "util.h"  // clear_memory
#include "str_manip.h"  // rev_comp
#include "hts.h"  // generic sequencer classes

using namespace Rcpp;







/*
 Sample for read lengths.

 If providing custom read lengths, make sure that the vector of their
 sampling weights is the same length as the vector of lengths.

 */
class PacBioReadLenSampler {

public:

    /* Initializers */
    PacBioReadLenSampler() {};
    // Using lognormal distribution:
    PacBioReadLenSampler(const double& scale_,
                         const double& sigma_,
                         const double& loc_,
                         const double& min_read_len_)
        : read_lens(),
          sampler(),
          distr(std::log(scale_), sigma_),
          use_distr(true),
          min_read_len(std::ceil(min_read_len_)),
          loc(loc_) {
        if (min_read_len < 1) min_read_len = 1;
    };
    // Using a vector of read lengths, each with a sampling probability:
    PacBioReadLenSampler(const std::vector<double>& read_probs_,
                         const std::vector<uint64>& read_lens_)
        : read_lens(read_lens_),
          sampler(read_probs_),
          distr(),
          use_distr(false),
          min_read_len(),
          loc() {
        if (read_probs_.size() != read_lens_.size()) {
            stop("Probability and read lengths vector should be the same length.");
        }
    };
    // Copy constructor
    PacBioReadLenSampler(const PacBioReadLenSampler& other)
        : read_lens(other.read_lens),
          sampler(other.sampler),
          distr(other.distr),
          use_distr(other.use_distr),
          min_read_len(other.min_read_len),
          loc(other.loc) {};
    // Assignment operator
    PacBioReadLenSampler& operator=(const PacBioReadLenSampler& other) {
        read_lens = other.read_lens;
        sampler = other.sampler;
        distr = other.distr;
        use_distr = other.use_distr;
        min_read_len = other.min_read_len;
        loc = other.loc;
        return *this;
    }

    uint64 sample(pcg64& eng);

private:

    std::vector<uint64> read_lens;      // optional vector of possible read lengths
    AliasSampler sampler;               // optional sampler that chooses from `read_lens`
    std::lognormal_distribution<double> distr; // optional if using a distribution
    bool use_distr;                     // Whether to sample using `distr` field
    double min_read_len;                // Minimum read length
    double loc;                         // Location parameter for lognormal distribution



};





/*
 Sample for number of passes.
 */
class PacBioPassSampler {

public:

    uint64 max_passes;
    std::vector<double> chi2_params_n;
    std::vector<double> chi2_params_s;


    /* Initializers */
    PacBioPassSampler() {};
    PacBioPassSampler(const uint64& max_passes_,
                      const std::vector<double>& chi2_params_n_,
                      const std::vector<double>& chi2_params_s_)
        : max_passes(max_passes_),
          chi2_params_n(chi2_params_n_),
          chi2_params_s(chi2_params_s_) {};
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


    void sample(uint64& split_pos,
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


        if ((static_cast<uint64>(wholes) & 1ULL) == 0ULL) {
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

};





class PacBioQualityError {

public:

    std::vector<double> sqrt_params;
    std::vector<double> norm_params;
    double prob_thresh;
    double prob_ins;
    double prob_del;
    double prob_subst;
    double min_exp;   // always initialize last

    /* Initializers */
    PacBioQualityError() {};
    PacBioQualityError(const std::vector<double>& sqrt_params_,
                       const std::vector<double>& norm_params_,
                       const double& prob_thresh_,
                       const double& prob_ins_,
                       const double& prob_del_,
                       const double& prob_subst_)
        : sqrt_params(sqrt_params_),
          norm_params(norm_params_),
          prob_thresh(prob_thresh_),
          prob_ins(prob_ins_),
          prob_del(prob_del_),
          prob_subst(prob_subst_),
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


    void sample(pcg64& eng,
                char& qual_left,
                char& qual_right,
                std::deque<uint64>& insertions,
                std::deque<uint64>& deletions,
                std::deque<uint64>& substitutions,
                const uint64& chrom_len,
                const uint64& read_length,
                const uint64& split_pos,
                const double& passes_left,
                const double& passes_right) {
        insertions.clear();
        deletions.clear();
        substitutions.clear();
        // Update sampling-error probabilities:
        update_probs(eng, passes_left, passes_right);
        // Update qualities
        fill_quals(qual_left, qual_right);
        // Now iterate through and update insertions, deletions, and substitutions:
        uint64 current_length = 0;
        uint64 chrom_pos = 0; // position on the read where events occur
        // Amount of extra (i.e., non-read) chromosome remaining:
        uint64 extra_space = chrom_len - read_length;
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
                    insertions.push_back(chrom_pos);
                    current_length++;
                    extra_space++;
                    if (current_length == split_pos) cum_probs = &cum_probs_right;
                }
                current_length++;
            } else if (u < cum_probs->at(1)) { // ----- deletion
                if (extra_space > 0) {
                    deletions.push_back(chrom_pos);
                    extra_space--;
                }
            } else { // ------------------------------- substitution
                substitutions.push_back(chrom_pos);
                current_length++;
            }
            chrom_pos++;
        }
        return;
    }



private:

    // Values that change:
    std::vector<double> cum_probs_left = std::vector<double>(3);
    std::vector<double> cum_probs_right = std::vector<double>(3);
    // Values that do not change:
    const uint64 max_qual = 93;
    const uint64 qual_start = static_cast<uint64>('!'); // quality of zero


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
         The first method is preferred over the second unless we're truncating
         all but the tail of the distribution:
         */
        if (lower_thresh < (norm_params[0] + 5 * norm_params[1])) {

            double p = R::pnorm5(a_bar, 0, 1, 1, 0);
            double u = runif_ab(eng, p, 1);

            double x = R::qnorm5(u, 0, 1, 1, 0);
            rnd = x * norm_params[1] + norm_params[0];

        } else {

            double u, x_bar, v;
            u = runif_01(eng);
            x_bar = std::sqrt(a_bar * a_bar  - 2 * std::log(1 - u));
            v = runif_01(eng);
            while (v > (x_bar / a_bar)) {
                u = runif_01(eng);
                x_bar = std::sqrt(a_bar * a_bar  - 2 * std::log(1 - u));
                v = runif_01(eng);
            }
            rnd = norm_params[1] * x_bar + norm_params[0];

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

    void update_probs(pcg64& eng,
                      const double& passes_left,
                      const double& passes_right);


    void fill_quals(char& qual_left, char& qual_right) {
        uint64 tmp_l = std::round(-10.0 * std::log10(cum_probs_left.back()));
        uint64 tmp_r = std::round(-10.0 * std::log10(cum_probs_right.back()));
        if (tmp_l > max_qual) tmp_l = max_qual;
        if (tmp_r > max_qual) tmp_r = max_qual;
        qual_left = static_cast<char>(tmp_l + qual_start);
        qual_right = static_cast<char>(tmp_r + qual_start);
        return;
    }





};







/*
 Template class to combine everything for PacBio sequencing of a single genome.
 (We will need multiple of these objects to chromosome a `HapSet` class.
 See `PacBioHaplotypes` class below.)

 `T` should be `HapGenome` or `RefGenome`

 */
template <typename T>
class PacBioOneGenome {
public:

    /* __ Samplers __ */
    // Samples read lengths:
    PacBioReadLenSampler len_sampler;
    // Samples numbers of passes over read:
    PacBioPassSampler pass_sampler;
    // Samples Pac Bio qualities and errors
    PacBioQualityError qe_sampler;


    /* __ Info __ */
    std::vector<uint64> chrom_reads;      // # reads per chromosome
    std::vector<uint64> chrom_lengths;    // genome-chromosome lengths
    const T* chromosomes;                 // pointer to `const T`
    std::string name;


    PacBioOneGenome() : chromosomes(nullptr) {};
    // Using lognormal distribution for read sizes:
    PacBioOneGenome(const T& chrom_object,
                    const double& scale_,
                    const double& sigma_,
                    const double& loc_,
                    const double& min_read_len_,
                    const uint64& max_passes_,
                    const std::vector<double>& chi2_params_n_,
                    const std::vector<double>& chi2_params_s_,
                    const std::vector<double>& sqrt_params_,
                    const std::vector<double>& norm_params_,
                    const double& prob_thresh_,
                    const double& prob_ins_,
                    const double& prob_del_,
                    const double& prob_subst_)
        : len_sampler(scale_, sigma_, loc_, min_read_len_),
          pass_sampler(max_passes_, chi2_params_n_, chi2_params_s_),
          qe_sampler(sqrt_params_, norm_params_, prob_thresh_, prob_ins_,
          prob_del_, prob_subst_),
          chrom_reads(),
          chrom_lengths(chrom_object.chrom_sizes()),
          chromosomes(&chrom_object),
          name(chrom_object.name) {};

    // Using vectors of read lengths and sampling weight for read lengths:
    PacBioOneGenome(const T& chrom_object,
                    const std::vector<double>& read_probs_,
                    const std::vector<uint64>& read_lens_,
                    const uint64& max_passes_,
                    const std::vector<double>& chi2_params_n_,
                    const std::vector<double>& chi2_params_s_,
                    const std::vector<double>& sqrt_params_,
                    const std::vector<double>& norm_params_,
                    const double& prob_thresh_,
                    const double& prob_ins_,
                    const double& prob_del_,
                    const double& prob_subst_)
        : len_sampler(read_probs_, read_lens_),
          pass_sampler(max_passes_, chi2_params_n_, chi2_params_s_),
          qe_sampler(sqrt_params_, norm_params_, prob_thresh_, prob_ins_,
          prob_del_, prob_subst_),
          chrom_reads(),
          chrom_lengths(chrom_object.chrom_sizes()),
          chromosomes(&chrom_object),
          name(chrom_object.name) {};

    PacBioOneGenome(const PacBioOneGenome& other)
        : len_sampler(other.len_sampler),
          pass_sampler(other.pass_sampler),
          qe_sampler(other.qe_sampler),
          chrom_reads(other.chrom_reads),
          chrom_lengths(other.chrom_lengths),
          chromosomes(other.chromosomes),
          name(other.name) {};



    void add_n_reads(const uint64& n_reads) {

        std::vector<double> probs_(chrom_lengths.begin(), chrom_lengths.end());
        chrom_reads = reads_per_group(n_reads, probs_);

        return;
    }

    // Add one read string (with 4 lines: ID, chromosome, "+", quality) to a FASTQ pool
    // `U` should be a std::string or std::vector<char>
    template <typename U>
    void one_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng);
    // Overloaded for when we input a haplotype chromosome stored as string
    template <typename U>
    void one_read(const std::string& chrom, const uint64& chrom_i,
                  std::vector<U>& fastq_pools, pcg64& eng);
    /*
     Same as above, but for a duplicate. It's assumed that `one_read` has been
     run once before.
     */
    template <typename U>
    void re_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng);
    template <typename U>
    void re_read(const std::string& chrom, const uint64& chrom_i,
                 std::vector<U>& fastq_pools, pcg64& eng);


private:

    uint64 split_pos = 0;
    double passes_left = 0;
    double passes_right = 0;
    char qual_left = '!';
    char qual_right = '!';
    uint64 read_chrom_space = 1;
    std::string read = std::string(1000, 'N');
    // Maps nucleotide char to integer from 0 to 3
    std::vector<uint8> nt_map = sequencer::nt_map;
    /*
    Maps nucleotide char integer (i.e., output from nt_map) to string of chars
    to sample from for a mismatch
    */
    std::vector<std::string> mm_nucleos = sequencer::mm_nucleos;
    // To store indel locations, where each vector will be of length 2 if paired==true
    std::deque<uint64> insertions = std::deque<uint64>(0);
    std::deque<uint64> deletions = std::deque<uint64>(0);
    std::deque<uint64> substitutions = std::deque<uint64>(0);
    uint64 chrom_ind = 0;
    uint64 read_length = 0;
    uint64 read_start = 0;

    // Append quality and read to fastq pool
    template <typename U>
    void append_pool(U& fastq_pool, pcg64& eng);
    template <typename U>
    void append_pool(const std::string& chrom, U& fastq_pool, pcg64& eng);


};


typedef PacBioOneGenome<RefGenome> PacBioReference;
typedef PacBioOneGenome<HapGenome> PacBioOneHaplotype;




/*
 To process a `HapSet` object, I need to wrap PacBioOneHaplotype inside
 another class.
 */
class PacBioHaplotypes {

public:

    const HapSet* haplotypes;                         // pointer to `const HapSet`
    std::vector<std::vector<uint64>> n_reads_vc;    // # reads per haplotype and chromosome
    std::vector<PacBioOneHaplotype> read_makers;      // makes PacBio reads
    std::vector<double> hap_probs;                  // probs of sampling haplotypes

    /* Initializers */
    PacBioHaplotypes() : haplotypes(nullptr) {};
    // Using lognormal distribution for read sizes:
    PacBioHaplotypes(const HapSet& hap_set,
                   const std::vector<double>& haplotype_probs,
                   const double& scale_,
                   const double& sigma_,
                   const double& loc_,
                   const double& min_read_len_,
                   const uint64& max_passes_,
                   const std::vector<double>& chi2_params_n_,
                   const std::vector<double>& chi2_params_s_,
                   const std::vector<double>& sqrt_params_,
                   const std::vector<double>& norm_params_,
                   const double& prob_thresh_,
                   const double& prob_ins_,
                   const double& prob_del_,
                   const double& prob_subst_)
        : haplotypes(&hap_set),
          n_reads_vc(),
          read_makers(),
          hap_probs(haplotype_probs),
          hap(0),
          chr(0),
          hap_chrom_seq() {

        /*
         Fill `read_makers` field:
         */
        uint64 n_haps = haplotypes->size();
        read_makers.reserve(n_haps);
        for (uint64 i = 0; i < n_haps; i++) {
            read_makers.push_back(
                PacBioOneHaplotype(hap_set[i],
                                 scale_, sigma_, loc_, min_read_len_,
                                 max_passes_, chi2_params_n_, chi2_params_s_,
                                 sqrt_params_, norm_params_, prob_thresh_,
                                 prob_ins_, prob_del_, prob_subst_)
            );
        }

    };
    // Using vectors of read lengths and sampling weight for read lengths:
    PacBioHaplotypes(const HapSet& hap_set,
                   const std::vector<double>& haplotype_probs,
                   const std::vector<double>& read_probs_,
                   const std::vector<uint64>& read_lens_,
                   const uint64& max_passes_,
                   const std::vector<double>& chi2_params_n_,
                   const std::vector<double>& chi2_params_s_,
                   const std::vector<double>& sqrt_params_,
                   const std::vector<double>& norm_params_,
                   const double& prob_thresh_,
                   const double& prob_ins_,
                   const double& prob_del_,
                   const double& prob_subst_)
        : haplotypes(&hap_set),
          n_reads_vc(),
          read_makers(),
          hap_probs(haplotype_probs),
          hap(0),
          chr(0),
          hap_chrom_seq() {

        /*
         Fill `read_makers` field:
         */
        uint64 n_haps = haplotypes->size();
        read_makers.reserve(n_haps);
        for (uint64 i = 0; i < n_haps; i++) {
            read_makers.push_back(
                PacBioOneHaplotype(hap_set[i],
                                 read_probs_, read_lens_,
                                 max_passes_, chi2_params_n_, chi2_params_s_,
                                 sqrt_params_, norm_params_, prob_thresh_,
                                 prob_ins_, prob_del_, prob_subst_)
            );
        }

    };

    // Copy constructor
    PacBioHaplotypes(const PacBioHaplotypes& other)
        : haplotypes(other.haplotypes),
          n_reads_vc(other.n_reads_vc),
          read_makers(other.read_makers),
          hap_probs(other.hap_probs),
          hap(other.hap),
          chr(other.chr),
          hap_chrom_seq(other.hap_chrom_seq) {};


    // Add info on # reads
    void add_n_reads(const uint64& n_reads) {

        uint64 n_haps = haplotypes->size();

        // split # reads by haplotype
        std::vector<uint64> hap_reads = reads_per_group(n_reads, hap_probs);
        // splitting by chromosome, too:
        for (uint64 v = 0; v < n_haps; v++) {
            std::vector<double> chrom_probs;
            for (const HapChrom& vc : (*haplotypes)[v].chromosomes) {
                chrom_probs.push_back(vc.size());
            }
            n_reads_vc.push_back(reads_per_group(hap_reads[v], chrom_probs));
        }

        // Fill `read_makers` field:
        for (uint64 i = 0; i < n_haps; i++) {
            read_makers[i].add_n_reads(hap_reads[i]);
        }

        return;
    }



    // `one_read` method
    template <typename U>
    void one_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng);
    // `re_read` method (for duplicates)
    template <typename U>
    void re_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng);



private:

    // Haplotype to create read from.
    uint64 hap;
    // Chromosome to create read from.
    uint64 chr;
    // String for haplotype chromosome. It's saved to make things faster.
    std::string hap_chrom_seq;


};







#endif
