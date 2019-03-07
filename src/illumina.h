#ifndef __GEMINO_ILLUMINA_H
#define __GEMINO_ILLUMINA_H



#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // distributions

#include "gemino_types.h"  // uint32
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // runif_01
#include "alias_sampler.h"  // AliasSampler
#include "util.h"  // clear_memory
#include "str_manip.h"  // rev_comp
#include "sequencer.h"  // other classes

using namespace Rcpp;








/*
 Sample for quality score when they vary by nucleotide.
 You'll want one of these objects for each position on the read.
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
 It then adds that error and outputs the quality string.
 For paired reads, you'll need two of these objects.
 */
class IlluminaQualityError {

public:

    std::vector<IllQualPos> by_nt; // One IllQualPos object per nucleotide type
    uint32 read_length;

    IlluminaQualityError() {};
    IlluminaQualityError(const std::vector<std::vector<std::vector<double>>>& probs_,
                         const std::vector<std::vector<std::vector<uint8>>>& quals_,
                         const double& ins_prob_,
                         const double& del_prob_)
        : by_nt(), read_length(probs_[0].size()),
          qual_prob_map(), nt_map(256, 4U), mm_nucleos(4),
          ins_prob(ins_prob_), del_prob(del_prob_),
          qual_start(static_cast<uint8>('!')) {

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

    // Sample for insertion and deletion positions
    void fill_indels(std::deque<uint32>& insertions,
                     std::deque<uint32>& deletions,
                     pcg64& eng) {
        uint32 i = 0, r = 0;
        double u;
        insertions.clear();
        deletions.clear();
        while (r < read_length) {
            u = runif_01(eng);
            if (u > (ins_prob + del_prob)) {
                r++;
            } else if (u > ins_prob) {
                deletions.push_back(i);
            } else {
                insertions.push_back(i);
                r += 2;
            }
            i++;
        }
        return;
    }

    // Fill read and quality strings
    void fill_read_qual(std::string& read,
                        std::string& qual,
                        std::deque<uint32>& insertions,
                        std::deque<uint32>& deletions,
                        pcg64& eng) const {
        uint32 seq_space = read_length + insertions.size() - deletions.size();
        if (read.size() != seq_space) stop("read.size() != seq_space");
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
    // Indel probabilities
    double ins_prob;
    double del_prob;
    /*
     Starting value of qualities (for use in converting integers to chars
     (e.g., 0 to '!'))
     */
    uint8 qual_start;

};



struct IlluminaFragLengths {

    std::gamma_distribution<double> distr;
    uint32 min;
    uint32 max;

    IlluminaFragLengths() {}
    IlluminaFragLengths(const double& shape, const double& scale,
                        const uint32& min_, const uint32& max_)
        : distr(shape, scale), min(min_), max(max_) {}

    uint32 sample(pcg64& eng) {
        double len_ = std::round(distr(eng));
        uint32 len = static_cast<uint32>(len_);
        if (len < min) len = min;
        if (len > max) len = max;
        return len;
    }

};



/*
 Incomplete template class to combine everything for Illumina sequencing.

 `T` should be `VarGenome` or `RefGenome`

 For full classes based on this template, you need to add a way to get read lengths.
 For Illumina, it will be a constant read length, and for long-read sequencing, it
 will be simply the fragment length.
 */
template <typename T>
class Illumina_t {
public:

    // __Sampler__                         __What it returns__
    TableSampler seqs;                  // index for which genome-sequence to sequence
    IlluminaFragLengths frag_lengths;   // fragment lengths
    IlluminaQualityError qual_errors1;  // Illumina qualities and errors, read 1
    IlluminaQualityError qual_errors2;  // Illumina qualities and errors, read 2


    // Info:
    std::vector<uint32> seq_lengths;    // genome-sequence lengths
    const T* sequences;                 // pointer to `const T`
    uint32 read_length;                 // Length of reads
    bool paired;                        // Boolean for whether to do paired-end reads

    Illumina_t() {};
    // For paired-end reads:
    Illumina_t(const T& seq_object,
          const double& frag_len_shape,
          const double& frag_len_scale,
          const uint32& frag_len_min,
          const uint32& frag_len_max,
          const std::vector<std::vector<std::vector<double>>>& mis_probs1,
          const std::vector<std::vector<std::vector<uint8>>>& quals1,
          const double& ins_prob1,
          const double& del_prob1,
          const std::vector<std::vector<std::vector<double>>>& mis_probs2,
          const std::vector<std::vector<std::vector<uint8>>>& quals2,
          const double& ins_prob2,
          const double& del_prob2)
        : seqs(),
          frag_lengths(frag_len_shape, frag_len_scale, frag_len_min, frag_len_max),
          qual_errors1(mis_probs1, quals1, ins_prob1, del_prob1),
          qual_errors2(mis_probs2, quals2, ins_prob2, del_prob2),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(mis_probs1[0].size()),
          paired(true) {
        if (mis_probs1[0].size() != mis_probs2[0].size()) {
            stop("In Illumina_t constr., read lengths for R1 and R2 don't match.");
        }
        this->construct_seqs();
    }
    // Single-end reads
    Illumina_t(const T& seq_object,
               const double& frag_len_shape,
               const double& frag_len_scale,
               const uint32& frag_len_min,
               const uint32& frag_len_max,
               const std::vector<std::vector<std::vector<double>>>& mis_probs,
               const std::vector<std::vector<std::vector<uint8>>>& quals,
               const double& ins_prob,
               const double& del_prob)
        : seqs(),
          frag_lengths(frag_len_shape, frag_len_scale, frag_len_min, frag_len_max),
          qual_errors1(mis_probs, quals, ins_prob, del_prob),
          qual_errors2(),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          read_length(mis_probs[0].size()),
          paired(false) {
        this->construct_seqs();
    }


private:

    std::deque<uint32> insertions;
    std::deque<uint32> deletions;

    void construct_seqs() {
        std::vector<double> probs_;
        probs_.reserve(seq_lengths.size());
        for (uint i = 0; i < seq_lengths.size(); i++) {
            probs_.push_back(static_cast<double>(seq_lengths[i]));
        }
        seqs = TableSampler(probs_);
    }


    /*

     LEFT OFF:  TURN THE BELOW FUNCTION INTO ONE THAT SAMPLES WHICH SEQUENCE TO USE,
     SAMPLE FOR INDELS, THEN SAMPLE WHERE TO PUT THE FRAGMENT
     [making sure the fragment size isn't less than the length required for the read:
     read_length + sum(indel sizes)]

     */

    // Sample a starting location for a fragment
    inline uint32 frag_loc_sample(const uint32& frag_len,
                                  const uint32& seq_len,
                                  pcg64& eng) {
        if (frag_len >= seq_len) return 0U;
        double u = runif_01(eng);
        uint32 pos = static_cast<uint32>(u * (seq_len - frag_len + 1));
        return pos;

        /*
         seq_ind = this->seqs.sample(eng);
         uint32 seq_len = (*(this->sequences))[seq_ind].size();

         frag_len = this->frag_lengths.sample(eng);
         if (frag_len > seq_len) frag_len = seq_len;
         frag_start = this->frag_loc_sample(frag_len, seq_len, eng);
         */
    }


};





#endif
