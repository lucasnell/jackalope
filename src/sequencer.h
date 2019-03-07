#ifndef __GEMINO_SEQUENCER_H
#define __GEMINO_SEQUENCER_H



#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng
#include <string>  // string class
#include <random>  // distributions

#include "gemino_types.h"  // uint32
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "pcg.h"  // runif_01
#include "table_sampler.h"  // TableSampler
#include "util.h"  // clear_memory
#include "str_manip.h"  // rev_comp

using namespace Rcpp;





/*
 Way to store info for the sequence identifier line in FASTQ file.
 It's organized as such:
 @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
    <read>:<is filtered>:<control number>:<sample number>
 (where the newline and tab are ' ')
 (from https://help.basespace.illumina.com/articles/descriptive/fastq-files/)
 Example: "@SIM:1:FCX:1:15:6329:1045 1:N:0:2"
 */
struct SequenceIdentifierInfo {

    // (Characters allowed: a–z, A–Z, 0–9 and underscore) Instrument ID:
    std::string instrument;
    // (Numerical) Run number on instrument:
    uint32 run_number;
    // (Characters allowed: a–z, A–Z, 0–9):
    std::string flowcell_ID;
    // (Numerical) Lane number:
    uint32 lane;
    // (Numerical) Tile number:
    uint32 tile;
    // (Numerical) Run number on instrument:
    uint32 x_pos;
    // (Numerical) X coordinate of cluster:
    uint32 y_pos;
    // (Numerical) Read number. 1 can be single read or Read 2 of paired-end:
    uint32 read;
    // (Y or N) Y if the read is filtered (did not pass), N otherwise:
    std::string is_filtered;
    // (Numerical) 0 when none of the control bits are on, otherwise even number:
    uint32 control_number;
    // (Numerical) Sample number from sample sheet:
    uint32 sample_number;

    SequenceIdentifierInfo() {}
    SequenceIdentifierInfo(
        const std::string& instrument_,
        const uint32& run_number_,
        const std::string& flowcell_ID_,
        const uint32& lane_,
        const uint32& tile_,
        const uint32& x_pos_,
        const uint32& y_pos_,
        const uint32& read_,
        const std::string& is_filtered_,
        const uint32& control_number_,
        const uint32& sample_number_
    )
        : instrument(instrument_), run_number(run_number_),
          flowcell_ID(flowcell_ID_), lane(lane_),
          tile(tile_), x_pos(x_pos_),
          y_pos(y_pos_), read(read_),
          is_filtered(is_filtered_), control_number(control_number_),
          sample_number(sample_number_) {};

    std::string get_id_line() {
        std::string out = "@";
        out += instrument + ':' + std::to_string(run_number) + ':' + flowcell_ID + ':';
        out += std::to_string(lane) + ':' + std::to_string(tile) + ':';
        out += std::to_string(x_pos) + ':' + std::to_string(y_pos) + ' ';
        out += std::to_string(read) + ':' + is_filtered + ':';
        out += std::to_string(control_number) + ':' + std::to_string(sample_number);
        return out;
    }

};








/*
 Sample for fragment lengths.
 This class is used to avoid having to store a probability for every single
 potential length, if the lengths can extend to something like 1 Mbp.

 It does weighted sampling using a vector of probabilities for each bin of lengths.
 Once it chooses a bin, it then uniformly samples within that bin if bin length > 1.
 */
class FragLenSampler {

public:

    uint32 region_len;                  // length of each region
    TableSampler sampler;               // sampler that chooses a region

    FragLenSampler() {};
    FragLenSampler(const std::vector<double>& probs,
                    const uint32& region_len_)
        : region_len(region_len_), sampler(probs) {}
    // Copy constructor
    FragLenSampler(const FragLenSampler& other)
        : region_len(other.region_len), sampler(other.sampler) {};
    // Assignment operator
    FragLenSampler& operator=(const FragLenSampler& other) {
        region_len = other.region_len;
        sampler = other.sampler;
        return *this;
    }

    uint32 sample(pcg64& eng) const {
        uint32 rnd = sampler.sample(eng);
        uint32 len_ = region_len * rnd;
        if (region_len > 1) len_ += static_cast<uint32>(runif_01(eng) * region_len);
        return len_;
    }

};





/*
 Sample for quality score when they vary by position on read.
 It does sampling using a normal distribution.
 The Mason sampler had separate quality rates for mismatched and non-mismatched
 read nucleotides.
 */
class QualitySampler {

public:

    uint32 region_len;                  // length of each region
    std::vector<double> means;          // means for each region
    std::vector<double> sds;            // SDs for each region

    QualitySampler() : region_len(), means(), sds(), distr(0, 1) {};
    QualitySampler(const std::vector<double>& means_,
                   const std::vector<double>& sds_,
                   const uint32& region_len_)
        : region_len(region_len_), means(means_), sds(sds_), distr(0, 1) {}

    QualitySampler& operator=(const QualitySampler& other) {
        region_len = other.region_len;
        means = other.means;
        sds = other.sds;
        return *this;
    }



    char sample(uint32 pos, pcg64& eng) const {
        // rnd ~ N(0, 1)
        double rnd = distr(eng);
        // Adjust pos to point to inside `means` and `sds` vectors
        pos /= region_len;
        if (pos >= means.size()) {
            Rcout << "`pos` is past QualitySampler bounds; using last value" << std::endl;
            pos = means.size() - 1;
        }
        // Adjust so rnd ~ N(means[pos], sds[pos])
        rnd *= sds[pos];
        rnd += means[pos];
        // Limit quality to [0, 41] (Illumina 1.8+ allows up to 41 instead of 40)
        if (rnd < 0) rnd = 0;
        if (rnd > 41) rnd = 41;
        // Convert to int, then char
        char q = static_cast<char>('!' + static_cast<int>(rnd));
        return q;
    }

private:

    mutable std::normal_distribution<double> distr;

};






/*
 Sample for sequencing error when rates vary by position on read.
 You should check that mismatch, insertion, and deletion probs never sum to > 1 outside
 this function!
 */
class SequencingErrorSampler {

public:

    uint32 region_len;                 // length of each region
    /*
     For each region, cumulative sum of Pr(match), Pr(mismatch), Pr(insertion)
     [deletion assumed otherwise]:
     */
    std::vector<std::vector<double>> probs;

    SequencingErrorSampler() : region_len(), probs() {};
    SequencingErrorSampler(const std::vector<double>& mis_probs,
                           const std::vector<double>& ins_probs,
                           const std::vector<double>& del_probs,
                           const uint32& region_len_)
        : region_len(region_len_), probs() {
        if (mis_probs.size() != ins_probs.size() ||
            mis_probs.size() != del_probs.size()) {
            stop("all input probs must be the same length.");
        }

        probs.reserve(mis_probs.size());
        std::vector<double> new_probs(3);

        for (uint32 i = 0; i < mis_probs.size(); i++) {
            double mis_i = mis_probs[i];
            double ins_i = ins_probs[i];
            double del_i = del_probs[i];
            double mat_i = 1 - (mis_i + ins_i + del_i);
            new_probs[0] = mat_i;
            new_probs[1] = mat_i + mis_i;
            new_probs[2] = mat_i + mis_i + ins_i;
            probs.push_back(new_probs);
        }

    }

    SequencingErrorSampler& operator=(const SequencingErrorSampler& other) {
        region_len = other.region_len;
        probs = other.probs;
        return *this;
    }



    /*
     Sample for error events for a read.
     Returns 0 for match, 1 for mismatch, 2 for insertion, 3 for deletion.
     Only 8-bit integer used to save on memory.
     */
    inline uint8 sample(uint32 pos, pcg64& eng) const {

        // Adjust pos to point to inside `probs` vector
        pos /= region_len;
        if (pos >= probs.size()) {
            Rcout << "`pos` is past SequencingErrorSampler bounds; using last value";
            Rcout << std::endl;
            pos = probs.size() - 1;
        }

        uint8 out;
        double u = runif_01(eng);
        if (u < probs[pos][0]) { // match
            out = 0;
        } else if (u < probs[pos][1]) { // mismatch
            out = 1;
        } else if (u < probs[pos][2]) { // insertion
            out = 2;
        } else {  // deletion
            out = 3;
        }

        return out;
    }

};



/*
 Incomplete template class to combine everything for whole-genome sequencing.

 `T` should be `VarGenome` or `RefGenome`

 For full classes based on this template, you need to add a way to get read lengths.
 For Illumina, it will be a constant read length, and for long-read sequencing, it
 will be simply the fragment length.
 */
template <typename T>
class WGS_t {
public:

    // __Sampler__                         __What it returns__
    TableSampler seqs;                  // index for which genome-sequence to sequence
    FragLenSampler frag_lengths;        // fragment lengths
    SequencingErrorSampler errors;      // sequencing errors
    TableSampler ins_lengths;           // insertion lengths
    TableSampler del_lengths;           // deletion lengths
    QualitySampler quals;               // qualities for matches
    QualitySampler mis_quals;           // qualities for mismatches

    // Info:
    std::vector<uint32> seq_lengths;    // genome-sequence lengths
    const T* sequences;                 // pointer to `const T`
    std::vector<std::string> mm;        // mismatch vector

    WGS_t()
        : seqs(), frag_lengths(), errors(), ins_lengths(), del_lengths(), quals(),
          mis_quals(), seq_lengths(), sequences(), mm(256, ""), bases("TCAGN") {};
    WGS_t(const T& seq_object,
          const std::vector<double>& frag_len_probs,
          const uint32& frag_len_region_len,
          const std::vector<double>& mis_probs,
          const std::vector<double>& ins_probs,
          const std::vector<double>& del_probs,
          const uint32& error_region_len,
          const std::vector<double>& ins_length_probs,
          const std::vector<double>& del_length_probs,
          const std::vector<double>& qual_means,
          const std::vector<double>& qual_sds,
          const uint32& qual_region_len,
          const std::vector<double>& mis_qual_means,
          const std::vector<double>& mis_qual_sds,
          const uint32& mis_qual_region_len)
        : seqs(),
          frag_lengths(frag_len_probs, frag_len_region_len),
          errors(mis_probs, ins_probs, del_probs, error_region_len),
          ins_lengths(ins_length_probs),
          del_lengths(del_length_probs),
          quals(qual_means, qual_sds, qual_region_len),
          mis_quals(mis_qual_means, mis_qual_sds, mis_qual_region_len),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          mm(256, ""),
          bases("TCAGN") {

        std::vector<double> probs_;
        probs_.reserve(seq_lengths.size());
        for (uint i = 0; i < seq_lengths.size(); i++) {
            probs_.push_back(static_cast<double>(seq_lengths[i]));
        }
        seqs = TableSampler(probs_);

        // For mismatches. Converts input char to uint, then uses that as index.
        mm[static_cast<uint32>('T')] = "CAG";
        mm[static_cast<uint32>('C')] = "TAG";
        mm[static_cast<uint32>('A')] = "TCG";
        mm[static_cast<uint32>('G')] = "TCA";
        mm[static_cast<uint32>('N')] = "NNN";
    }


protected:

    std::string bases;

    // Sample a starting location for a fragment
    inline uint32 frag_loc_sample(const uint32& frag_len,
                                  const uint32& seq_len,
                                  pcg64& eng) {
        if (frag_len >= seq_len) return 0U;
        double u = runif_01(eng);
        uint32 pos = static_cast<uint32>(u * (seq_len - frag_len + 1));
        return pos;
    }

    // Adjust one character for a mismatch
    void mismatch(char& input, pcg64& eng) {
        const std::string& str(mm[static_cast<uint32>(input)]);
        uint32 ind = static_cast<uint32>(runif_01(eng) * str.size());
        input = str[ind];
        return;
    }

    // Adjust string and qualities for an insertion
    void insertion(const sint32& size, const uint32& read_pos,
                   std::string& read, std::string& qual,
                   pcg64& eng) {
        uint32 size_ = static_cast<uint32>(size);
        std::string nts;
        nts.reserve(size_);
        for (uint32 j = 0; j < size_; j++) {
            uint32 k = runif_01(eng) * 4;
            nts += this->bases[k];
            qual += this->mis_quals.sample(read_pos, eng);
        }
        read.insert(read_pos + 1, nts);
        return;
    }


};




/*
 Sequencing object for whole genome sequencing of one reference/variant using Illumina
 */

template <typename T>
class IlluminaWGS_t: public WGS_t<T> {
public:

    uint32 read_length;
    uint32 n_reads;

    // Constructors:
    IlluminaWGS_t() : WGS_t<T>(), read_length(), n_reads() {}
    IlluminaWGS_t(const T& seq_object,
                  const std::vector<double>& frag_len_probs,
                  const uint32& frag_len_region_len,
                  const std::vector<double>& mis_probs,
                  const std::vector<double>& ins_probs,
                  const std::vector<double>& del_probs,
                  const uint32& error_region_len,
                  const std::vector<double>& ins_length_probs,
                  const std::vector<double>& del_length_probs,
                  const std::vector<double>& qual_means,
                  const std::vector<double>& qual_sds,
                  const uint32& qual_region_len,
                  const std::vector<double>& mis_qual_means,
                  const std::vector<double>& mis_qual_sds,
                  const uint32& mis_qual_region_len,
                  const uint32& read_length_,
                  const bool& paired) :
        WGS_t<T>(seq_object, frag_len_probs, frag_len_region_len,
                 mis_probs, ins_probs, del_probs, error_region_len,
                 ins_length_probs, del_length_probs,
                 qual_means, qual_sds, qual_region_len,
                 mis_qual_means, mis_qual_sds, mis_qual_region_len),
        read_length(read_length_), n_reads(paired ? uint32(2) : uint32(1)) {};
    // Copy constructor
    IlluminaWGS_t<T>(const IlluminaWGS_t<T>& other)
        : WGS_t<T>(other),
          read_length(other.read_length),
          n_reads(other.n_reads) {};
    // Assignment operator
    IlluminaWGS_t<T>& operator=(const IlluminaWGS_t<T>& other) {
        WGS_t<T>::operator=(other);
        this->read_length = other.read_length;
        this->n_reads = other.n_reads;
        return *this;
    }


    // Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality)
    std::vector<std::string> one_read(pcg64& eng, SequenceIdentifierInfo& ID_info) {

        // Initiate objects
        uint32 seq_ind;
        uint32 frag_len;
        uint32 frag_start;
        std::vector<std::string> reads;
        std::vector<std::string> quals;
        std::vector<std::vector<uint8>> err_codes;
        std::vector<std::vector<sint32>> err_lengths;
        std::vector<uint32> read_seq_spaces;
        /*
         `read_seq_spaces` above is the "sequence space": read length + sum(indel lengths)
         */

        /*
        Fill all objects but `reads` and `quals` (those have vector-lengths set, but
        their inner strings are not set):
        */
        fill_info(seq_ind, frag_len, frag_start, reads, quals, err_codes,
                  err_lengths, read_seq_spaces, eng);

        std::vector<std::string> output(this->n_reads);

        // Boolean for whether we take the reverse side first:
        bool reverse = runif_01(eng) < 0.5;
        for (uint32 i = 0; i < this->n_reads; i++) {
            ID_info.read = i + 1;
            std::string& read(reads[i]);
            std::string& qual(quals[i]);
            qual.reserve(std::max(read_seq_spaces[i], this->read_length));
            // Read starting location:
            uint32 start = frag_start;
            if (reverse) start += (frag_len - read_seq_spaces[i]);

            // Now fill `read` from `sequences` field:
            (*(this->sequences))[seq_ind].fill_seq(read, start, read_seq_spaces[i]);

            // Reverse-complement `read` if taking reverse side:
            if (reverse) rev_comp(read);

            // Adjust for errors:
            adjust_for_errors(err_codes[i], err_lengths[i], read, qual, eng);

            // If doing paired reads, the second one should be the reverse of the first
            reverse = !reverse;

            // Combine into 4 lines of output per read:
            output[i] = ID_info.get_id_line() + '\n' + read + "\n+\n" + qual;
        }

        return output;

    }


protected:


    /*
     Fill in vectors of error codes and sizes, plus add up how much total
     "sequence space" the read takes up.
     The latter is the read length plus the sum of all indel-error lengths.
     This number is increased by deletions (bc they skip sequence) and decreased
     by insertions.
     */
    void fill_errors(std::vector<uint8>& err_codes,
                     std::vector<sint32>& err_lengths,
                     uint32& read_seq_space,
                     const uint32& read_len,
                     const uint32& frag_len,
                     pcg64& eng) {

        // Restart all these:
        if (err_codes.size() > 0) err_codes.clear();
        if (err_lengths.size() > 0) err_lengths.clear();
        read_seq_space = read_len;

        // Reserve memory (`* 1.1` in case of more deletions than insertions)
        err_codes.reserve(static_cast<uint32>(read_len * 1.1));
        err_lengths.reserve(static_cast<uint32>(read_len * 1.1));

        uint32 i = 0;
        while (i < read_len) {
            uint8 err_code = this->errors.sample(i, eng);
            err_codes.push_back(err_code);
            if (err_code == 0) { // ----------------------------------- match
                i++;
                err_lengths.push_back(0);
            } else if (err_code == 1) { // ---------------------------- mismatch
                i++;
                err_lengths.push_back(0);
            } else if (err_code == 2) { // ---------------------------- insertion
                sint32 size_ = static_cast<sint32>(this->ins_lengths.sample(eng) + 1);
                i += (size_ + 1);
                // Adjust for an insert going beyond read length:
                if (i > read_len) size_ -= (i - read_len);
                err_lengths.push_back(size_);
                read_seq_space -= size_;
            } else if (err_code == 3) { // ---------------------------- deletion
                sint32 size_ = static_cast<sint32>(this->del_lengths.sample(eng) + 1);
                size_ *= -1;
                err_lengths.push_back(size_);
                read_seq_space -= size_;
                /*
                If, because of this deletion, the read goes past the fragment length,
                then we set the sequence space to the fragment length and stop here.
                */
                if (read_seq_space >= frag_len) {
                    read_seq_space = frag_len;
                    break;
                }
            }
        }

        return;
    }



    /*
     Do most things except for actually filling in the sequence.
     Filling in the sequence will be dependent on the class used in the `sequences`
     field (i.e., as the `T` argument to this template)
     */
    void fill_info(uint32& seq_ind,
                   uint32& frag_len,
                   uint32& frag_start,
                   std::vector<std::string>& reads,
                   std::vector<std::string>& quals,
                   std::vector<std::vector<uint8>>& err_codes,
                   std::vector<std::vector<sint32>>& err_lengths,
                   std::vector<uint32>& read_seq_spaces,
                   pcg64& eng) {

        seq_ind = this->seqs.sample(eng);
        uint32 seq_len = (*(this->sequences))[seq_ind].size();

        frag_len = this->frag_lengths.sample(eng);
        if (frag_len > seq_len) frag_len = seq_len;
        frag_start = this->frag_loc_sample(frag_len, seq_len, eng);

        uint32 read_len = this->read_length;
        if (read_len > frag_len) read_len = frag_len;

        // Create read and quality strings
        reads = std::vector<std::string>(this->n_reads);
        quals = std::vector<std::string>(this->n_reads);
        /*
         Start and fill vectors of error codes and sizes, to know how much
         "sequence space" (read length + sum(indel lengths)) we'll need.
         */
        err_codes = std::vector<std::vector<uint8>>(this->n_reads);
        err_lengths = std::vector<std::vector<sint32>>(this->n_reads);
        read_seq_spaces = std::vector<uint32>(this->n_reads);
        for (uint32 i = 0; i < this->n_reads; i++) {
            this->fill_errors(err_codes[i], err_lengths[i], read_seq_spaces[i],
                              read_len, frag_len, eng);
        }

        return;

    }


    // Adjust for sequencing errors:
    void adjust_for_errors(const std::vector<uint8>& err_codes,
                           const std::vector<sint32>& err_lengths,
                           std::string& read,
                           std::string& qual,
                           pcg64& eng) {
        uint32 read_pos = 0;
        for (uint32 j = 0; j < err_codes.size(); j++) {
            const uint8& err_code(err_codes[j]);
            if (err_code == 0) { // -------------------- match
                qual += this->quals.sample(read_pos, eng);
                read_pos++;
            } else if (err_code == 1) { // ------------- mismatch
                this->mismatch(read[read_pos], eng);
                qual += this->mis_quals.sample(read_pos, eng);
                read_pos++;
            } else if (err_code == 2) { // ------------- insertion
                const sint32& size_(err_lengths[j]);
                // Adjust string and qualities:
                this->insertion(size_, read_pos, read, qual, eng);
                read_pos += (size_ + 1);
            } else if (err_code == 3) { // ------------- deletion
                const sint32& size_(err_lengths[j]);
                // (below, it won't try to delete past read size)
                read.erase(read_pos, size_);
            }
        }
    }
};


typedef IlluminaWGS_t<RefGenome> ReferenceIlluminaWGS;



/*
 To process a `VarSet` object, I need to wrap IlluminaWGS_t<VarGenome> inside
 another class.
 */
class VariantIlluminaWGS {

public:

    const VarSet* variants;                 // pointer to `const VarSet`
    TableSampler variant_sampler;           // chooses which variant to use
    IlluminaWGS_t<VarGenome> read_maker;    // makes Illumina reads

    VariantIlluminaWGS() {}
    VariantIlluminaWGS(const VarSet& var_set,
                       const std::vector<double>& variant_probs,
                       const std::vector<double>& frag_len_probs,
                       const uint32& frag_len_region_len,
                       const std::vector<double>& mis_probs,
                       const std::vector<double>& ins_probs,
                       const std::vector<double>& del_probs,
                       const uint32& error_region_len,
                       const std::vector<double>& ins_length_probs,
                       const std::vector<double>& del_length_probs,
                       const std::vector<double>& qual_means,
                       const std::vector<double>& qual_sds,
                       const uint32& qual_region_len,
                       const std::vector<double>& mis_qual_means,
                       const std::vector<double>& mis_qual_sds,
                       const uint32& mis_qual_region_len,
                       const uint32& read_length_,
                       const bool& paired) :
        variants(&var_set),
        variant_sampler(variant_probs),
        read_maker(var_set[0], frag_len_probs, frag_len_region_len,
                   mis_probs, ins_probs, del_probs, error_region_len,
                   ins_length_probs, del_length_probs,
                   qual_means, qual_sds, qual_region_len,
                   mis_qual_means, mis_qual_sds, mis_qual_region_len,
                   read_length_, paired) {};


    /*
     -------------
     `one_read` methods
     -------------
     */
    // If only providing rng and id info, sample for a variant, then make read(s):
    std::vector<std::string> one_read(pcg64& eng,
                                      SequenceIdentifierInfo& ID_info) {
        uint32 var = variant_sampler.sample(eng);
        read_maker.sequences = &((*variants)[var]);
        std::vector<std::string> out = read_maker.one_read(eng, ID_info);
        return out;
    }
    // If you provide a specific variant, then make read(s) from that:
    std::vector<std::string> one_read(const uint32& var,
                                      pcg64& eng,
                                      SequenceIdentifierInfo& ID_info) {
        read_maker.sequences = &((*variants)[var]);
        std::vector<std::string> out = read_maker.one_read(eng, ID_info);
        return out;
    }


};






/*
 Sequencing object for whole genome sequencing of one reference/variant using long-read
 sequencing (Nanopore/PacBio)
*/

template <typename T>
class LongReadWGS_t: public WGS_t<T> {
public:

    // Constructors:
    LongReadWGS_t() : WGS_t<T>() {}
    LongReadWGS_t(const T& seq_object,
                  const std::vector<double>& frag_len_probs,
                  const uint32& frag_len_region_len,
                  const std::vector<double>& mis_probs,
                  const std::vector<double>& ins_probs,
                  const std::vector<double>& del_probs,
                  const uint32& error_region_len,
                  const std::vector<double>& ins_length_probs,
                  const std::vector<double>& del_length_probs,
                  const std::vector<double>& qual_means,
                  const std::vector<double>& qual_sds,
                  const uint32& qual_region_len,
                  const std::vector<double>& mis_qual_means,
                  const std::vector<double>& mis_qual_sds,
                  const uint32& mis_qual_region_len) :
    WGS_t<T>(seq_object, frag_len_probs, frag_len_region_len,
             mis_probs, ins_probs, del_probs, error_region_len,
             ins_length_probs, del_length_probs,
             qual_means, qual_sds, qual_region_len,
             mis_qual_means, mis_qual_sds, mis_qual_region_len) {};
    // Copy constructor
    LongReadWGS_t<T>(const LongReadWGS_t<T>& other)
        : WGS_t<T>(other) {};
    // Assignment operator
    LongReadWGS_t<T>& operator=(const LongReadWGS_t<T>& other) {
        WGS_t<T>::operator=(other);
        return *this;
    }


    // Sample one read string (with 4 lines: ID, sequence, "+", quality)
    std::string one_read(pcg64& eng, SequenceIdentifierInfo& ID_info) {

        uint32 seq_ind = this->seqs.sample(eng);
        uint32 seq_len = (*(this->sequences))[seq_ind].size();

        uint32 frag_len = this->frag_lengths.sample(eng);
        if (frag_len > seq_len) frag_len = seq_len;
        uint32 frag_start = this->frag_loc_sample(frag_len, seq_len, eng);

        std::string read;
        std::string qual;
        std::string output;  // This will combine both, plus have ID and separator line
        ID_info.read = 1;

        qual.reserve(static_cast<uint32>(frag_len * 1.1)); // `*1.1` in case of insertions

        // Boolean for whether we take the reverse side:
        bool reverse = runif_01(eng) < 0.5;

        // Now fill `read` from `sequences` field:
        (*(this->sequences))[seq_ind].fill_seq(read, frag_start, frag_len);

        // Reverse-complement `read` if taking reverse side:
        if (reverse) rev_comp(read);

        // Add in sequencing errors:
        add_errors(read, qual, eng);

        // Combine into 4 lines of output per read:
        output = ID_info.get_id_line() + '\n' + read + "\n+\n" + qual;

        return output;

    }


protected:



    // Add sequencing errors to read:
    void add_errors(std::string& read, std::string& qual, pcg64& eng) {

        uint32 read_pos = 0;
        while (read_pos < read.size()) {
            uint8 err_code = this->errors.sample(read_pos, eng);
            if (err_code == 0) { // -------------------------------------- match
                qual += this->quals.sample(read_pos, eng);
                read_pos++;
            } else if (err_code == 1) { // ------------------------------- mismatch
                this->mismatch(read[read_pos], eng);
                qual += this->mis_quals.sample(read_pos, eng);
                read_pos++;
            } else if (err_code == 2) { // ------------------------------- insertion
                sint32 size_ = static_cast<sint32>(this->ins_lengths.sample(eng) + 1);
                // Adjust string and qualities:
                this->insertion(size_, read_pos, read, qual, eng);
                read_pos += (size_ + 1);
            } else if (err_code == 3) { // ------------------------------- deletion
                sint32 size_ = static_cast<sint32>(this->del_lengths.sample(eng) + 1);
                size_ *= -1;
                // (below, it won't try to delete past read size)
                read.erase(read_pos, size_);
            }
        }

        return;
    }
};

typedef LongReadWGS_t<RefGenome> ReferenceLongReadWGS;




/*
 To process a `VarSet` object, I need to wrap LongReadWGS_t<VarGenome> inside
 another class.
 */
class VariantLongReadWGS {

public:

    const VarSet* variants;                 // pointer to `const VarSet`
    TableSampler variant_sampler;           // chooses which variant to use
    LongReadWGS_t<VarGenome> read_maker;    // makes Illumina reads

    VariantLongReadWGS() {}
    VariantLongReadWGS(const VarSet& var_set,
                       const std::vector<double>& variant_probs,
                       const std::vector<double>& frag_len_probs,
                       const uint32& frag_len_region_len,
                       const std::vector<double>& mis_probs,
                       const std::vector<double>& ins_probs,
                       const std::vector<double>& del_probs,
                       const uint32& error_region_len,
                       const std::vector<double>& ins_length_probs,
                       const std::vector<double>& del_length_probs,
                       const std::vector<double>& qual_means,
                       const std::vector<double>& qual_sds,
                       const uint32& qual_region_len,
                       const std::vector<double>& mis_qual_means,
                       const std::vector<double>& mis_qual_sds,
                       const uint32& mis_qual_region_len) :
        variants(&var_set),
        variant_sampler(variant_probs),
        read_maker(var_set[0], frag_len_probs, frag_len_region_len,
                   mis_probs, ins_probs, del_probs, error_region_len,
                   ins_length_probs, del_length_probs,
                   qual_means, qual_sds, qual_region_len,
                   mis_qual_means, mis_qual_sds, mis_qual_region_len) {};

    /*
     -------------
     `one_read` methods
     -------------
     */
    // If only providing rng and id info, sample for a variant, then make read:
    std::string one_read(pcg64& eng, SequenceIdentifierInfo& ID_info) {
        uint32 var = variant_sampler.sample(eng);
        read_maker.sequences = &((*variants)[var]);
        std::string out = read_maker.one_read(eng, ID_info);
        return out;
    }
    // If you provide a specific variant, then make read from that:
    std::string one_read(const uint32& var, pcg64& eng, SequenceIdentifierInfo& ID_info) {
        read_maker.sequences = &((*variants)[var]);
        std::string out = read_maker.one_read(eng, ID_info);
        return out;
    }


};











#endif
