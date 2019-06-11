#ifndef __JACKAL_ILLUMINA_H
#define __JACKAL_ILLUMINA_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // distributions
#include <fstream> // for writing FASTQ files
#include "zlib.h"  // for writing to compressed FASTQ
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // runif_01
#include "alias_sampler.h"  // AliasSampler
#include "hts.h"  // generic sequencing class


using namespace Rcpp;




/*
 ----------------------------------------------------------------
 ----------------------------------------------------------------
 ----------------------------------------------------------------

 CREATING READS

 ----------------------------------------------------------------
 ----------------------------------------------------------------
 ----------------------------------------------------------------
 */


// Basic information to construct reads
struct IlluminaReadConstrInfo {
    uint64 read_length;
    uint64 seq_ind;
    uint64 frag_len;
    uint64 frag_start;
    std::vector<std::string> reads;
    std::vector<std::string> quals;
    std::vector<uint64> read_seq_spaces;
    std::string barcode;


    IlluminaReadConstrInfo() {}
    IlluminaReadConstrInfo(const bool& paired,
                           const uint64& read_length_,
                           const std::string barcode_)
        : read_length(read_length_),
          seq_ind(0),
          frag_len(0),
          frag_start(0),
          reads(),
          quals(),
          read_seq_spaces(),
          barcode(barcode_) {

        if (paired) {
            reads = std::vector<std::string>(2, std::string(read_length, 'N'));
            quals = std::vector<std::string>(2);
            read_seq_spaces = std::vector<uint64>(2);
        } else {
            reads = std::vector<std::string>(1, std::string(read_length, 'N'));
            quals = std::vector<std::string>(1);
            read_seq_spaces = std::vector<uint64>(1);
        }
    }
    IlluminaReadConstrInfo(const IlluminaReadConstrInfo& other)
        : read_length(other.read_length), seq_ind(other.seq_ind),
          frag_len(other.frag_len), frag_start(other.frag_start),
          reads(other.reads), quals(other.quals),
          read_seq_spaces(other.read_seq_spaces), barcode(other.barcode) {};
};




/*
 Sample for quality score when they vary by position on read.
 You'll want one of these objects for each nucleotide.
 */
class IllQualPos {

public:

    std::vector<AliasSampler> samplers;
    std::vector<std::vector<uint8>> quals;
    uint64 read_length;

    IllQualPos() {};
    IllQualPos(const std::vector<std::vector<double>>& probs_,
               const std::vector<std::vector<uint8>>& quals_)
        : samplers(), quals(quals_), read_length(quals_.size()) {

        if (probs_.size() != read_length) {
            stop("In IllQualPos construct, probs_.size() != quals_.size()");
        }

        samplers.reserve(read_length);
        quals.reserve(read_length);
        for (uint64 i = 0; i < read_length; i++) {
            samplers.push_back(AliasSampler(probs_[i]));
        }

    }
    IllQualPos(const IllQualPos& other)
        : samplers(other.samplers), quals(other.quals), read_length(other.read_length) {};

    IllQualPos& operator=(const IllQualPos& other) {
        samplers = other.samplers;
        quals = other.quals;
        read_length = other.read_length;
        return *this;
    }

    // Sample for a quality
    uint8 sample(const uint64& pos,
                 pcg64& eng) const {
        uint64 k = samplers[pos].sample(eng);
        return quals[pos][k];
    }



};




/*
 Sample for (1) quality score and (2) whether a mismatch error is produced.
 It does NOT sample for indels, but it does add indels to the read string.
 It also adds mismatch errors to the read string and outputs the quality string.
 For paired reads, you'll need two of these objects.
 */
class IlluminaQualityError {

public:

    std::vector<IllQualPos> by_nt; // One IllQualPos object per nucleotide type

    IlluminaQualityError() {};
    IlluminaQualityError(const std::vector<std::vector<std::vector<double>>>& probs_,
                         const std::vector<std::vector<std::vector<uint8>>>& quals_)
        : by_nt(),
          qual_prob_map() {

        uint64 read_length(probs_[0].size());

        if (probs_.size() != 4 || quals_.size() != 4) {
            stop("All probs and quals for IlluminaQualityError must be of length 4");
        }

        by_nt.reserve(4);
        // For making the vector to map qualities to probabilities of mismatches:
        uint8 max_qual = 0;
        for (uint64 i = 0; i < 4; i++) {
            if (probs_[i].size() != read_length) {
                stop("In IlluminaQualityError construct, all probs' lengths not equal");
            }
            if (quals_[i].size() != read_length) {
                stop("In IlluminaQualityError construct, all quals' lengths not equal");
            }
            by_nt.push_back(IllQualPos(probs_[i], quals_[i]));
            for (const std::vector<uint8>& qvec : quals_[i]) {
                uint8 max_ij = *std::max_element(qvec.begin(), qvec.end());
                if (max_ij > max_qual) max_qual = max_ij;
            }
        }

        qual_prob_map.reserve(max_qual+1);  // `+1` bc we're using qualities as indices
        qual_prob_map.push_back(1);
        for (uint64 q = 1; q < (static_cast<uint64>(max_qual)+1ULL); q++) {
            double prob = std::pow(10, static_cast<double>(q) / -10.0);
            qual_prob_map.push_back(prob);
        }

    }

    IlluminaQualityError(const IlluminaQualityError& other)
        : by_nt(other.by_nt),
          qual_prob_map(other.qual_prob_map) {};


    /*
     Fill read and quality strings.
     Because small fragments could cause the read length to be less than normal,
     I'm requiring it as input to this function.
     `read` should already be sized appropriately before this function.
     */
    void fill_read_qual(std::string& read,
                        std::string& qual,
                        std::deque<uint64>& insertions,
                        std::deque<uint64>& deletions,
                        pcg64& eng) const {

        double mis_prob, u;
        uint8 nt_ind, qint;
        /*
         Add indels:
         */
        uint64 seq_pos = read.size() - 1ULL;
        while (!insertions.empty() || !deletions.empty()) {
            if (!insertions.empty() && seq_pos == insertions.back()) {
                char c = alias_sampler::bases[static_cast<uint64>(runif_01(eng) * 4.0)];
                read.insert(seq_pos + 1, 1, c);
                insertions.pop_back();
            } else if (!deletions.empty() && seq_pos == deletions.back()) {
                read.erase(seq_pos, 1);
                deletions.pop_back();
            }
            if (seq_pos == 0) break;
            seq_pos--;
        }
        if (qual.size() != read.size()) qual.resize(read.size());
        /*
         Add mismatches:
         */
        for (uint64 pos = 0; pos < read.size(); pos++) {
            char& nt(read[pos]);
            nt_ind = nt_map[nt];
            /*
             For all values except for T, C, A, or G, it'll return random quality less
             than 10. This is what ART does.
             */
            if (nt_ind > 3) {
                qint = runif_01(eng) * 10 + qual_start;
                qual[pos] = static_cast<char>(qint);
                nt = 'N';
                continue;
            }
            /*
             Otherwise, qualities are based on the nucleotide and position,
             and Pr(mismatch) is proportional to quality:
             */
            qint = by_nt[nt_ind].sample(pos, eng);
            mis_prob = qual_prob_map[qint];
            qint += qual_start;
            qual[pos] = static_cast<char>(qint);
            u = runif_01(eng);
            if (u < mis_prob) {
                const std::string& mm_str(mm_nucleos[nt_ind]);
                nt = mm_str[static_cast<uint64>(runif_01(eng) * 3.0)];
            }
        }

        return;
    }

private:
    // Maps quality integer to probability of mismatch double:
    std::vector<double> qual_prob_map;
    // Maps nucleotide char to integer from 0 to 3
    std::vector<uint8> nt_map = sequencer::nt_map;
    /*
     Maps nucleotide char integer (i.e., output from nt_map) to string of chars
     to sample from for a mismatch
     */
    std::vector<std::string> mm_nucleos = sequencer::mm_nucleos;
    /*
     Starting value of qualities (for use in converting integers to chars
     (e.g., 0 to '!'))
     */
    uint8 qual_start = static_cast<uint8>('!');

};







/*
 Template class to combine everything for Illumina sequencing of a single genome.
 (We will need multiple of these objects to sequence a `VarSet` class.
  See `IlluminaVariants` class below.)

 `T` should be `VarGenome` or `RefGenome`

 */
template <typename T>
class IlluminaOneGenome {
public:

    /* __ Samplers __ */
    // Samples index for which genome-sequence to sequence
    AliasSampler seq_sampler;
    // Samples Illumina qualities and errors, one `IlluminaQualityError` for each read
    std::vector<IlluminaQualityError> qual_errors;
    // Samples fragment lengths:
    std::gamma_distribution<double> frag_lengths;   // fragment lengths


    /* __ Info __ */
    std::vector<uint64> seq_lengths;    // genome-sequence lengths
    const T* sequences;                 // pointer to `const T`
    uint64 read_length;                 // Length of reads
    bool paired;                        // Boolean for whether to do paired-end reads
    bool matepair;                      // Boolean for whether to do mate-pair reads
    std::vector<double> ins_probs;      // Per-base prob. of an insertion, reads 1 and 2
    std::vector<double> del_probs;      // Per-base prob. of a deletion, reads 1 and 2
    std::string name;

    IlluminaOneGenome() : sequences(nullptr) {};
    // For paired-end reads:
    IlluminaOneGenome(const T& seq_object,
                      const bool& matepair_,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint64& frag_len_min_,
                      const uint64& frag_len_max_,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                      const std::vector<std::vector<std::vector<uint8>>>& quals1,
                      const double& ins_prob1,
                      const double& del_prob1,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                      const std::vector<std::vector<std::vector<uint8>>>& quals2,
                      const double& ins_prob2,
                      const double& del_prob2,
                      const std::string& barcode)
        : seq_sampler(),
          qual_errors(),
          frag_lengths(frag_len_shape, frag_len_scale),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(qual_probs1[0].size()),
          paired(true),
          matepair(matepair_),
          ins_probs(2),
          del_probs(2),
          name(seq_object.name),
          insertions(2),
          deletions(2),
          frag_len_min(frag_len_min_),
          frag_len_max(frag_len_max_),
          constr_info(paired, read_length, barcode) {
              if (qual_probs1[0].size() != qual_probs2[0].size()) {
                  std::string err = "In IlluminaOneGenome constr., read lengths for ";
                  err += "R1 and R2 don't match.";
                  stop(err.c_str());
              }
              qual_errors = {IlluminaQualityError(qual_probs1, quals1),
                                   IlluminaQualityError(qual_probs2, quals2)};
              construct_seqs();
              ins_probs[0] = ins_prob1;
              ins_probs[1] = ins_prob2;
              del_probs[0] = del_prob1;
              del_probs[1] = del_prob2;
          };
    // Single-end reads
    IlluminaOneGenome(const T& seq_object,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint64& frag_len_min_,
                      const uint64& frag_len_max_,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs,
                      const std::vector<std::vector<std::vector<uint8>>>& quals,
                      const double& ins_prob,
                      const double& del_prob,
                      const std::string& barcode)
        : seq_sampler(),
          qual_errors{IlluminaQualityError(qual_probs, quals)},
          frag_lengths(frag_len_shape, frag_len_scale),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(qual_probs[0].size()),
          paired(false),
          matepair(false),
          ins_probs(1),
          del_probs(1),
          name(seq_object.name),
          insertions(1),
          deletions(1),
          frag_len_min(frag_len_min_),
          frag_len_max(frag_len_max_),
          constr_info(paired, read_length, barcode) {
              construct_seqs();
              ins_probs[0] = ins_prob;
              del_probs[0] = del_prob;
          };

    IlluminaOneGenome(const IlluminaOneGenome& other)
        : seq_sampler(other.seq_sampler),
          qual_errors(other.qual_errors),
          frag_lengths(other.frag_lengths),
          seq_lengths(other.seq_lengths),
          sequences(other.sequences),
          read_length(other.read_length),
          paired(other.paired),
          matepair(other.matepair),
          ins_probs(other.ins_probs),
          del_probs(other.del_probs),
          name(other.name),
          insertions(other.insertions),
          deletions(other.deletions),
          frag_len_min(other.frag_len_min),
          frag_len_max(other.frag_len_max),
          constr_info(other.constr_info) {};



    // Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality)
    // `U` should be a std::string or std::vector<char>
    template <typename U>
    void one_read(std::vector<U>& fastq_pools, pcg64& eng);

    /*
     Same as above, but for a duplicate. It's assumed that `one_read` has been
     run once before.
     */
    template <typename U>
    void re_read(std::vector<U>& fastq_pools, pcg64& eng);

    /*
     Add information about a RefGenome or VarGenome object
     This is used when making multiple samplers that share most info except for
     that related to the sequence object.
     */
    void add_seq_info(const T& seq_object, const std::string& barcode);




protected:

    // To store indel locations, where each vector will be of length 2 if paired==true
    std::vector<std::deque<uint64>> insertions;
    std::vector<std::deque<uint64>> deletions;
    // Bounds for fragment sizes:
    uint64 frag_len_min;
    uint64 frag_len_max;
    // Info to construct reads:
    IlluminaReadConstrInfo constr_info;


    // Construct sequence-sampling probabilities:
    void construct_seqs();


    // Sample for insertion and deletion positions
    void sample_indels(pcg64& eng);

    // Adjust sequence spaces
    void adjust_seq_spaces();


    /*
     Sample a sequence, indels, fragment length, and starting position for the fragment.
     Lastly, it sets the sequence spaces required for these reads.
     */
    void seq_indels_frag(pcg64& eng);


    /*
     Same as above, but for duplicates.
     This means skipping the sequence and fragment info parts.
     */
    void just_indels(pcg64& eng);



    /*
     Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality),
     then append that to the `fastq_pools` vector.
     This function does NOT do anything with fragments.
     That should be done outside this function.
     */
    template <typename U>
    void append_pools(std::vector<U>& fastq_pools, pcg64& eng);


};


typedef IlluminaOneGenome<RefGenome> IlluminaReference;
typedef IlluminaOneGenome<VarGenome> IlluminaOneVariant;



/*
 To process a `VarSet` object, I need to wrap IlluminaOneVariant inside
 another class.
 */
class IlluminaVariants {

public:

    const VarSet* variants;                         // pointer to `const VarSet`
    AliasSampler variant_sampler;                   // chooses which variant to use
    std::vector<IlluminaOneVariant> read_makers;    // makes Illumina reads
    bool paired;                                    // Boolean for paired-end reads

    IlluminaVariants() : variants(nullptr) {}

    /* Initializers */

    // For paired-end reads:
    IlluminaVariants(const VarSet& var_set,
                     const std::vector<double>& variant_probs,
                     const bool& matepair_,
                     const double& frag_len_shape,
                     const double& frag_len_scale,
                     const uint64& frag_len_min_,
                     const uint64& frag_len_max_,
                     const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                     const std::vector<std::vector<std::vector<uint8>>>& quals1,
                     const double& ins_prob1,
                     const double& del_prob1,
                     const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                     const std::vector<std::vector<std::vector<uint8>>>& quals2,
                     const double& ins_prob2,
                     const double& del_prob2,
                     std::vector<std::string> barcodes)
        : variants(&var_set),
          variant_sampler(variant_probs),
          read_makers(),
          paired(true),
          var(0) {

        if (barcodes.size() < var_set.size()) barcodes.resize(var_set.size(), "");

        /*
         Fill `read_makers` field:
         */
        uint64 n_vars = var_set.size();
        // Read maker for the first variant:
        IlluminaOneVariant read_maker1(var_set[0], matepair_,
                                       frag_len_shape, frag_len_scale,
                                       frag_len_min_, frag_len_max_,
                                       qual_probs1, quals1, ins_prob1, del_prob1,
                                       qual_probs2, quals2, ins_prob2, del_prob2,
                                       barcodes[0]);
        read_makers.reserve(n_vars);
        read_makers.push_back(read_maker1);
        for (uint64 i = 1; i < n_vars; i++) {
            read_makers.push_back(read_maker1);
            read_makers[i].add_seq_info(var_set[i], barcodes[i]);
        }

    };

    // Single-end reads
    IlluminaVariants(const VarSet& var_set,
                     const std::vector<double>& variant_probs,
                     const double& frag_len_shape,
                     const double& frag_len_scale,
                     const uint64& frag_len_min_,
                     const uint64& frag_len_max_,
                     const std::vector<std::vector<std::vector<double>>>& qual_probs,
                     const std::vector<std::vector<std::vector<uint8>>>& quals,
                     const double& ins_prob,
                     const double& del_prob,
                     std::vector<std::string> barcodes)
        : variants(&var_set),
          variant_sampler(variant_probs),
          read_makers(),
          paired(false),
          var(0) {

        if (barcodes.size() < var_set.size()) barcodes.resize(var_set.size(), "");

        /*
         Fill `read_makers` field:
         */
        uint64 n_vars = var_set.size();
        // Read maker for the first variant:
        IlluminaOneVariant read_maker1(var_set[0],
                                       frag_len_shape, frag_len_scale,
                                       frag_len_min_, frag_len_max_,
                                       qual_probs, quals, ins_prob, del_prob,
                                       barcodes[0]);
        read_makers.reserve(n_vars);
        read_makers.push_back(read_maker1);
        for (uint64 i = 1; i < n_vars; i++) {
            read_makers.push_back(read_maker1);
            read_makers[i].add_seq_info(var_set[i], barcodes[i]);
        }

    };

    IlluminaVariants(const IlluminaVariants& other)
        : variants(other.variants), variant_sampler(other.variant_sampler),
          read_makers(other.read_makers), paired(other.paired), var(other.var) {};


    /*
     -------------
     `one_read` methods
     -------------
     */
    // If only providing rng and id info, sample for a variant, then make read(s):
    template <typename U>
    void one_read(std::vector<U>& fastq_pools, pcg64& eng) {
        var = variant_sampler.sample(eng);
        read_makers[var].one_read<U>(fastq_pools, eng);
        return;
    }
    /*
     -------------
     `re_read` methods (for duplicates)
     -------------
     */
    template <typename U>
    void re_read(std::vector<U>& fastq_pools, pcg64& eng) {
        read_makers[var].re_read<U>(fastq_pools, eng);
        return;
    }



private:

    // Variant to sample from. It's saved in this class in case of duplicates.
    uint64 var;

};





#endif
