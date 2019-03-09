#ifndef __GEMINO_ILLUMINA_H
#define __GEMINO_ILLUMINA_H



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // distributions

#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // runif_01
#include "alias_sampler.h"  // AliasSampler
#include "table_sampler.h"  // TableSampler
#include "util.h"  // clear_memory
#include "str_manip.h"  // rev_comp
#include "sequencer.h"  // SequenceIdentifierInfo class

using namespace Rcpp;






// Basic information to construct reads
struct ReadConstructInfo {
    uint32 read_length;
    uint32 seq_ind;
    uint32 frag_len;
    uint32 frag_start;
    std::vector<std::string> reads;
    std::vector<std::string> quals;
    std::vector<uint32> read_seq_spaces;

    ReadConstructInfo() {}
    ReadConstructInfo(const bool& paired, const uint32& read_length_) {
        this->read_length = read_length_;
        if (paired) {
            reads = std::vector<std::string>(2);
            quals = std::vector<std::string>(2);
            read_seq_spaces = std::vector<uint32>(2);
        } else {
            reads = std::vector<std::string>(1);
            quals = std::vector<std::string>(1);
            read_seq_spaces = std::vector<uint32>(1);
        }
    }
};




/*
 Sample for quality score when they vary by position on read.
 You'll want one of these objects for each nucleotide.
 */
class IllQualPos {

public:

    std::vector<AliasSampler> samplers;
    std::vector<std::vector<uint8>> quals;
    uint32 read_length;

    IllQualPos() {};
    IllQualPos(const std::vector<std::vector<double>>& probs_,
               const std::vector<std::vector<uint8>>& quals_)
        : samplers(), quals(quals_), read_length(quals_.size()) {

        if (probs_.size() != read_length) {
            stop("In IllQualPos construct, probs_.size() != quals_.size()");
        }

        samplers.reserve(read_length);
        quals.reserve(read_length);
        for (uint32 i = 0; i < read_length; i++) {
            samplers.push_back(AliasSampler(probs_[i]));
        }

    }

    IllQualPos& operator=(const IllQualPos& other) {
        this->samplers = other.samplers;
        this->quals = other.quals;
        this->read_length = other.read_length;
        return *this;
    }

    // Sample for a quality
    uint8 sample(const uint32& pos,
                 pcg64& eng) const {
        uint32 k = samplers[pos].sample(eng);
        return quals[pos][k];
    }

    /*
     Sample for a quality and update the probability of a mismatch based on that quality
     */
    uint8 sample(const uint32& pos,
                 double& mis_prob,
                 const std::vector<double>& qual_prob_map,
                 pcg64& eng) const {

        uint32 k = samplers[pos].sample(eng);
        uint8 qual = quals[pos][k];
        mis_prob = qual_prob_map[qual];
        return qual;
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
          qual_prob_map(), nt_map(256, 4U), mm_nucleos(4),
          qual_start(static_cast<uint8>('!')) {

        uint32 read_length(probs_[0].size());

        if (probs_.size() != 4 || quals_.size() != 4) {
            stop("All probs and quals for IlluminaQualityError must be of length 4");
        }

        by_nt.reserve(4);
        // For making the vector to map qualities to probabilities of mismatches:
        uint8 max_qual = 0;
        for (uint32 i = 0; i < 4; i++) {
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
        for (uint32 q = 1; q < (max_qual+1); q++) {
            double prob = std::pow(10, static_cast<double>(q) / -10.0);
            qual_prob_map.push_back(prob);
        }

        // Filling map to go from nucleotide char to integer from 0 to 3
        // (returns 4 if not T, C, A, or G)
        for (uint32 i = 0; i < 4; i++) {
            nt_map[alias_sampler::bases[i]] = i;
        }
        mm_nucleos[nt_map[static_cast<uint8>('T')]] = "CAG";
        mm_nucleos[nt_map[static_cast<uint8>('C')]] = "TAG";
        mm_nucleos[nt_map[static_cast<uint8>('A')]] = "TCG";
        mm_nucleos[nt_map[static_cast<uint8>('G')]] = "TCA";

    }


    /*
     Fill read and quality strings.
     Because small fragments could cause the read length to be less than normal,
     I'm requiring it as input to this function.
     */
    void fill_read_qual(std::string& read,
                        std::string& qual,
                        std::deque<uint32>& insertions,
                        std::deque<uint32>& deletions,
                        const uint32& read_length,
                        pcg64& eng) const {

        uint32 seq_space = read_length + insertions.size() - deletions.size();
        if (qual.size() != read_length) qual.resize(read_length);
        double mis_prob, u;
        uint8 nt_ind, qint;
        /*
         Add indels:
         */
        uint32 seq_pos = seq_space - 1;
        while (seq_pos > 0) {
            if (!insertions.empty() && seq_pos == insertions.back()) {
                char c = alias_sampler::bases[static_cast<uint32>(runif_01(eng) * 4UL)];
                read.insert(seq_pos + 1, 1, c);
                insertions.pop_back();
            } else if (!deletions.empty() && seq_pos == deletions.back()) {
                read.erase(seq_pos, 1);
                deletions.pop_back();
            }
            seq_pos--;
        }
        if (read.size() != read_length) stop("read.size() != read_length after indels.");
        /*
         Add mismatches:
         */
        for (uint32 pos = 0; pos < read_length; pos++) {
            char& nt(read[pos]);
            nt_ind = nt_map[nt];
            /*
             For all values except for T, C, A, or G, it'll return random quality less
             than 10. This is what ART does.
             */
            if (nt_ind > 3) {
                qint = runif_01(eng) * 10 + this->qual_start;
                qual[pos] = static_cast<char>(qint);
                nt = 'N';
                continue;
            }
            /*
             Otherwise, qualities are based on the nucleotide and position,
             and Pr(mismatch) is proportional to quality:
             */
            qint = by_nt[nt_ind].sample(pos, mis_prob, this->qual_prob_map, eng) +
                this->qual_start;
            qual[pos] = static_cast<char>(qint);
            u = runif_01(eng);
            if (u < mis_prob) {
                const std::string& mm_str(mm_nucleos[nt_map[nt]]);
                nt = mm_str[runif_aabb(eng, 0U, 2U)];
            }
        }
    }

private:
    // Maps quality integer to probability of mismatch double:
    std::vector<double> qual_prob_map;
    // Maps nucleotide char to integer from 0 to 3
    std::vector<uint8> nt_map;
    /*
     Maps nucleotide char integer (i.e., output from nt_map) to string of chars
     to sample from for a mismatch
     */
    std::vector<std::string> mm_nucleos;
    /*
     Starting value of qualities (for use in converting integers to chars
     (e.g., 0 to '!'))
     */
    uint8 qual_start;

};







/*
 Template class to combine everything for Illumina sequencing.

 `T` should be `VarGenome` or `RefGenome`

 */
template <typename T>
class Illumina_t {
public:

    /* __ Samplers __ */
    // Samples index for which genome-sequence to sequence
    AliasSampler seqs;
    // Samples Illumina qualities and errors, one `IlluminaQualityError` for each read
    std::vector<IlluminaQualityError> qual_errors;
    // Samples fragment lengths:
    std::gamma_distribution<double> frag_lengths;   // fragment lengths


    /* __ Info __ */
    std::vector<uint32> seq_lengths;    // genome-sequence lengths
    const T* sequences;                 // pointer to `const T`
    uint32 read_length;                 // Length of reads
    bool paired;                        // Boolean for whether to do paired-end reads
    std::vector<double> ins_probs;      // Per-base prob. of an insertion, reads 1 and 2
    std::vector<double> del_probs;      // Per-base prob. of a deletion, reads 1 and 2

    Illumina_t() {};
    // For paired-end reads:
    Illumina_t(const T& seq_object,
               const double& frag_len_shape,
               const double& frag_len_scale,
               const uint32& frag_len_min_,
               const uint32& frag_len_max_,
               const std::vector<std::vector<std::vector<double>>>& qual_probs1,
               const std::vector<std::vector<std::vector<uint8>>>& quals1,
               const double& ins_prob1,
               const double& del_prob1,
               const std::vector<std::vector<std::vector<double>>>& qual_probs2,
               const std::vector<std::vector<std::vector<uint8>>>& quals2,
               const double& ins_prob2,
               const double& del_prob2)
        : seqs(),
          qual_errors(),
          frag_lengths(frag_len_shape, frag_len_scale),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(qual_probs1[0].size()),
          paired(true),
          ins_probs{ins_prob1, ins_prob2},
          del_probs{del_prob1, del_prob2},
          insertions(2),
          deletions(2),
          frag_len_min(frag_len_min_),
          frag_len_max(frag_len_max_),
          constr_info(paired, read_length) {
              if (qual_probs1[0].size() != qual_probs2[0].size()) {
                  stop("In Illumina_t constr., read lengths for R1 and R2 don't match.");
              }
              this->qual_errors = {IlluminaQualityError(qual_probs1, quals1),
                                   IlluminaQualityError(qual_probs2, quals2)};
              this->construct_seqs();
          };
    // Single-end reads
    Illumina_t(const T& seq_object,
               const double& frag_len_shape,
               const double& frag_len_scale,
               const uint32& frag_len_min_,
               const uint32& frag_len_max_,
               const std::vector<std::vector<std::vector<double>>>& qual_probs,
               const std::vector<std::vector<std::vector<uint8>>>& quals,
               const double& ins_prob,
               const double& del_prob)
        : seqs(),
          qual_errors{IlluminaQualityError(qual_probs, quals)},
          frag_lengths(frag_len_shape, frag_len_scale),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(qual_probs[0].size()),
          paired(false),
          ins_probs{ins_prob},
          del_probs{del_prob},
          insertions(1),
          deletions(1),
          frag_len_min(frag_len_min_),
          frag_len_max(frag_len_max_),
          constr_info(paired, read_length) {
              this->construct_seqs();
          };


    // Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality)
    void one_read(std::vector<std::string>& read_quals,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info) {

        /*
         Sample fragment info, and set the sequence space(s) required for these read(s).
         */
        this->seq_indels_frag(eng);

        // Fill the reads and qualities
        this->fill_strings(read_quals, eng, ID_info);

        return;
    }

    /*
     Same as above, but for a PCR duplicate. It's assumed that `one_read` has been
     run once before.
     */
    void re_read(std::vector<std::string>& read_quals,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info) {

        // Here I'm just re-doing indels bc it's a PCR duplicate.
        this->just_indels(eng);

        // Fill the reads and qualities
        this->fill_strings(read_quals, eng, ID_info);

        return;
    }

    /*
     Add information about a RefGenome or VarGenome object
     This is used when making multiple samplers that share most info except for
     that related to the sequence object.
     */
    void add_seq_info(const T& seq_object) {

        this->seq_lengths = seq_object.seq_sizes();
        this->sequences = &seq_object;

        std::vector<double> probs_;
        probs_.reserve(this->seq_lengths.size());
        for (uint i = 0; i < seq_lengths.size(); i++) {
            probs_.push_back(static_cast<double>(seq_lengths[i]));
        }
        this->seqs = AliasSampler(probs_);

        return;
    }


protected:

    // To store indel locations, where each vector will be of length 2 if paired==true
    std::vector<std::deque<uint32>> insertions;
    std::vector<std::deque<uint32>> deletions;
    // Bounds for fragment sizes:
    uint32 frag_len_min;
    uint32 frag_len_max;
    // Info to construct reads:
    ReadConstructInfo constr_info;


    // Construct sequence-sampling probabilities:
    inline void construct_seqs() {
        std::vector<double> probs_;
        probs_.reserve(this->seq_lengths.size());
        for (uint i = 0; i < this->seq_lengths.size(); i++) {
            probs_.push_back(static_cast<double>(this->seq_lengths[i]));
        }
        this->seqs = AliasSampler(probs_);
        return;
    }


    // Sample for insertion and deletion positions
    void sample_indels(pcg64& eng, const uint32& frag_len) {

        for (uint32 r = 0; r < insertions.size(); r++) {
            uint32 i = 0, read_size = 0;
            double u;
            std::deque<uint32>& ins(this->insertions[r]);
            std::deque<uint32>& del(this->deletions[r]);
            const double& ins_prob(ins_probs[r]);
            const double& del_prob(del_probs[r]);
            ins.clear();
            del.clear();
            while (read_size < read_length && read_size < frag_len) {
                u = runif_01(eng);
                if (u > (ins_prob + del_prob)) {
                    read_size++;
                } else if (u > ins_prob) {
                    del.push_back(i);
                } else {
                    ins.push_back(i);
                    read_size += 2;
                }
                i++;
            }
        }

        return;
    }


    /*
     Sample a sequence, indels, fragment length, and starting position for the fragment.
     Lastly, it sets the sequence spaces required for these reads.
     */
    inline void seq_indels_frag(pcg64& eng) {

        uint32& seq_ind(this->constr_info.seq_ind);
        uint32& frag_len(this->constr_info.frag_len);
        uint32& frag_start(this->constr_info.frag_start);

        // Sample sequence:
        seq_ind = this->seqs.sample(eng);
        uint32 seq_len = (*(this->sequences))[seq_ind].size();

        // Sample fragment length:
        frag_len = static_cast<uint32>(this->frag_lengths(eng));
        if (frag_len < frag_len_min) frag_len = frag_len_min;
        if (frag_len > frag_len_max) frag_len = frag_len_max;

        // Sample fragment starting position:
        if (frag_len >= seq_len) {
            frag_len = seq_len;
            frag_start = 0;
        } else {
            double u = runif_01(eng);
            frag_start = static_cast<uint32>(u * (seq_len - frag_len + 1));
        }
        // Sample indels:
        this->sample_indels(eng, frag_len);

        // Set sequence spaces:
        std::vector<uint32>& read_seq_spaces(this->constr_info.read_seq_spaces);
        uint32 len_ = std::min(read_length, frag_len);
        for (uint i = 0; i < read_seq_spaces.size(); i++) {
            read_seq_spaces[i] = len_ + insertions[i].size() -
                deletions[i].size();
        }

        return;
    }


    /*
     Same as above, but for PCR duplicates.
     This means skipping the sequence and fragment info parts.
     */
    inline void just_indels(pcg64& eng) {

        const uint32& frag_len(this->constr_info.frag_len);
        // Sample indels:
        this->sample_indels(eng, frag_len);

        // Set sequence spaces:
        std::vector<uint32>& read_seq_spaces(this->constr_info.read_seq_spaces);
        uint32 len_ = std::min(read_length, frag_len);
        for (uint i = 0; i < read_seq_spaces.size(); i++) {
            read_seq_spaces[i] = len_ + insertions[i].size() -
                deletions[i].size();
        }

        return;
    }



    /*
     Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality).
     This function does NOT do anything with fragments.
     That should be done outside this function.
     */
    void fill_strings(std::vector<std::string>& read_quals,
                      pcg64& eng,
                      SequenceIdentifierInfo& ID_info) {

        uint32 n_reads = this->ins_probs.size();
        if (read_quals.size() != n_reads) read_quals.resize(n_reads);

        // If we sampled a very small fragment, it'll reduce read length:
        uint32 real_read_length = std::min(read_length, this->constr_info.frag_len);

        // Boolean for whether we take the reverse side first:
        bool reverse = runif_01(eng) < 0.5;
        for (uint32 i = 0; i < n_reads; i++) {
            ID_info.read = i + 1;
            std::string& read(this->constr_info.reads[i]);
            std::string& qual(this->constr_info.quals[i]);

            // Read starting location:
            uint32 start = this->constr_info.frag_start;
            if (reverse) start += (this->constr_info.frag_len -
                this->constr_info.read_seq_spaces[i]);

            // Now fill `read` from `sequences` field:
            (*(this->sequences))[this->constr_info.seq_ind].fill_seq(read, start,
             this->constr_info.read_seq_spaces[i]);

            // Reverse-complement `read` if taking reverse side:
            if (reverse) rev_comp(read);

            // Sample mapping quality and add errors to read:
            qual_errors[i].fill_read_qual(read, qual, insertions[i], deletions[i],
                                          real_read_length, eng);

            // If doing paired reads, the second one should be the reverse of the first
            reverse = !reverse;

            // Combine into 4 lines of output per read:
            read_quals[i] = ID_info.get_id_line() + '\n' + read + "\n+\n" + qual;
        }

        return;
    }


};


typedef Illumina_t<RefGenome> ReferenceIllumina;



/*
 To process a `VarSet` object, I need to wrap Illumina_t<VarGenome> inside
 another class.
 */
class VariantIllumina {

public:

    const VarSet* variants;                         // pointer to `const VarSet`
    TableSampler variant_sampler;                   // chooses which variant to use
    std::vector<Illumina_t<VarGenome>> read_makers; // makes Illumina reads

    VariantIllumina() {}

    /*
     ------------------------
     Initializers

     Initialization doesn't include VarSet information, so it just starts `read_makers`
     with one item that includes all info but the `VarSet` info that will come later.
     ------------------------
     */

    // For paired-end reads:
    VariantIllumina(const VarSet& var_set,
                    const std::vector<double>& variant_probs,
                    const double& frag_len_shape,
                    const double& frag_len_scale,
                    const uint32& frag_len_min_,
                    const uint32& frag_len_max_,
                    const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                    const std::vector<std::vector<std::vector<uint8>>>& quals1,
                    const double& ins_prob1,
                    const double& del_prob1,
                    const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                    const std::vector<std::vector<std::vector<uint8>>>& quals2,
                    const double& ins_prob2,
                    const double& del_prob2)
        : variants(&var_set),
          variant_sampler(variant_probs),
          read_makers(),
          var(0) {

        /*
         Fill `read_makers` field:
         */
        uint32 n_vars = var_set.size();
        // Read maker for the first variant:
        Illumina_t<VarGenome> read_maker1(var_set[0],
                                          frag_len_shape, frag_len_scale,
                                          frag_len_min_, frag_len_max_,
                                          qual_probs1, quals1, ins_prob1, del_prob1,
                                          qual_probs2, quals2, ins_prob2, del_prob2);
        read_makers.reserve(n_vars);
        read_makers.push_back(read_maker1);
        for (uint32 i = 1; i < n_vars; i++) {
            read_makers.push_back(read_maker1);
            read_makers[i].add_seq_info(var_set[i]);
        }

    };

    // Single-end reads
    VariantIllumina(const VarSet& var_set,
                    const std::vector<double>& variant_probs,
                    const double& frag_len_shape,
                    const double& frag_len_scale,
                    const uint32& frag_len_min_,
                    const uint32& frag_len_max_,
                    const std::vector<std::vector<std::vector<double>>>& qual_probs,
                    const std::vector<std::vector<std::vector<uint8>>>& quals,
                    const double& ins_prob,
                    const double& del_prob)
        : variants(&var_set),
          variant_sampler(variant_probs),
          read_makers(),
          var(0) {

        /*
         Fill `read_makers` field:
         */
        uint32 n_vars = var_set.size();
        // Read maker for the first variant:
        Illumina_t<VarGenome> read_maker1(var_set[0],
                                          frag_len_shape, frag_len_scale,
                                          frag_len_min_, frag_len_max_,
                                          qual_probs, quals, ins_prob, del_prob);
        read_makers.reserve(n_vars);
        read_makers.push_back(read_maker1);
        for (uint32 i = 1; i < n_vars; i++) {
            read_makers.push_back(read_maker1);
            read_makers[i].add_seq_info(var_set[i]);
        }

    };


    /*
     -------------
     `one_read` methods
     -------------
     */
    // If only providing rng and id info, sample for a variant, then make read(s):
    void one_read(std::vector<std::string>& read_quals,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info) {
        this->var = variant_sampler.sample(eng);
        read_makers[this->var].one_read(read_quals, eng, ID_info);
        return;
    }
    // If you provide a specific variant, then make read(s) from that:
    void one_read(std::vector<std::string>& read_quals,
                  const uint32& var_,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info) {
        this->var = var_;
        read_makers[this->var].one_read(read_quals, eng, ID_info);
        return;
    }
    /*
     -------------
     `re_read` methods (for PCR duplicates)
     -------------
     */
    void re_read(std::vector<std::string>& read_quals,
                  pcg64& eng,
                  SequenceIdentifierInfo& ID_info) {
        read_makers[this->var].re_read(read_quals, eng, ID_info);
        return;
    }


private:

    // Variant to sample from. It's saved in this class in case of PCR duplicates.
    uint32 var;

};







#endif
