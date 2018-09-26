#ifndef __GEMINO_SEQUENCER_H
#define __GEMINO_SEQUENCER_H


/*
 ==========================================================================
                         Mason - A Read Simulator
 ==========================================================================
 Copyright (c) 2006-2016, Knut Reinert, FU Berlin
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
     * Neither the name of Knut Reinert or the FU Berlin nor the names of
       its contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.

 ==========================================================================
 Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ==========================================================================

 Edited for use in gemino by Lucas Nell, 2018

 */

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


namespace sequencer {
    std::string bases = "TCAGN";
}



class SeqOptions {

public:

    double read_length;

    // Sequencing-error arguments
    double prob_insert;
    double prob_delete;
    double prob_mismatch_scale;
    double prob_mismatch;
    double prob_mismatch_begin;
    double prob_mismatch_end;
    double position_raise;

    // Quality arguments
    double mean_qual_begin;
    double mean_qual_end;
    double sd_qual_begin;
    double sd_qual_end;
    double mean_mismatch_qual_begin;
    double mean_mismatch_qual_end;
    double sd_mismatch_qual_begin;
    double sd_mismatch_qual_end;

    // Constructor
    SeqOptions() {

        read_length = 100;

        prob_insert = 0.001;
        prob_delete = 0.001;
        prob_mismatch_scale = 1.0;
        prob_mismatch = 0.004;
        prob_mismatch_begin = 0.002;
        prob_mismatch_end = 0.012;
        position_raise = 0.66;

        mean_qual_begin = 40;
        mean_qual_end = 39.5;
        sd_qual_begin = 0.05;
        sd_qual_end = 10;
        mean_mismatch_qual_begin = 39.5;
        mean_mismatch_qual_end = 30;
        sd_mismatch_qual_begin = 3;
        sd_mismatch_qual_end = 15;

    }

    SeqOptions(const List& opt_list) {

        SeqOptions defaults;

        if (opt_list.containsElementNamed("read_length") ) {
            double rl_check = static_cast<double>(opt_list["read_length"]);
            if (rl_check < 1) stop("read length cannot be < 1");
            read_length = opt_list["read_length"];
        } else read_length = defaults.read_length;

        if (opt_list.containsElementNamed("prob_insert") ) {
            prob_insert = opt_list["prob_insert"];
        } else prob_insert = defaults.prob_insert;

        if (opt_list.containsElementNamed("prob_delete") ) {
            prob_delete = opt_list["prob_delete"];
        } else prob_delete = defaults.prob_delete;

        if (opt_list.containsElementNamed("prob_mismatch_scale") ) {
            prob_mismatch_scale = opt_list["prob_mismatch_scale"];
        } else prob_mismatch_scale = defaults.prob_mismatch_scale;

        if (opt_list.containsElementNamed("prob_mismatch") ) {
            prob_mismatch = opt_list["prob_mismatch"];
        } else prob_mismatch = defaults.prob_mismatch;

        if (opt_list.containsElementNamed("prob_mismatch_begin") ) {
            prob_mismatch_begin = opt_list["prob_mismatch_begin"];
        } else prob_mismatch_begin = defaults.prob_mismatch_begin;

        if (opt_list.containsElementNamed("prob_mismatch_end") ) {
            prob_mismatch_end = opt_list["prob_mismatch_end"];
        } else prob_mismatch_end = defaults.prob_mismatch_end;

        if (opt_list.containsElementNamed("position_raise") ) {
            position_raise = opt_list["position_raise"];
        } else position_raise = defaults.position_raise;

        if (opt_list.containsElementNamed("mean_qual_begin") ) {
            mean_qual_begin = opt_list["mean_qual_begin"];
        } else mean_qual_begin = defaults.mean_qual_begin;

        if (opt_list.containsElementNamed("mean_qual_end") ) {
            mean_qual_end = opt_list["mean_qual_end"];
        } else mean_qual_end = defaults.mean_qual_end;

        if (opt_list.containsElementNamed("sd_qual_begin") ) {
            sd_qual_begin = opt_list["sd_qual_begin"];
        } else sd_qual_begin = defaults.sd_qual_begin;

        if (opt_list.containsElementNamed("sd_qual_end") ) {
            sd_qual_end = opt_list["sd_qual_end"];
        } else sd_qual_end = defaults.sd_qual_end;

        if (opt_list.containsElementNamed("mean_mismatch_qual_begin") ) {
            mean_mismatch_qual_begin = opt_list["mean_mismatch_qual_begin"];
        } else mean_mismatch_qual_begin = defaults.mean_mismatch_qual_begin;

        if (opt_list.containsElementNamed("mean_mismatch_qual_end") ) {
            mean_mismatch_qual_end = opt_list["mean_mismatch_qual_end"];
        } else mean_mismatch_qual_end = defaults.mean_mismatch_qual_end;

        if (opt_list.containsElementNamed("sd_mismatch_qual_begin") ) {
            sd_mismatch_qual_begin = opt_list["sd_mismatch_qual_begin"];
        } else sd_mismatch_qual_begin = defaults.sd_mismatch_qual_begin;

        if (opt_list.containsElementNamed("sd_mismatch_qual_end") ) {
            sd_mismatch_qual_end = opt_list["sd_mismatch_qual_end"];
        } else sd_mismatch_qual_end = defaults.sd_mismatch_qual_end;

    }

};


class SeqValues {

public:

    uint32 read_length;
    std::vector<double> mismatch_probs;
    double prob_insert;
    double prob_delete;
    std::vector<double> qual_mean;
    std::vector<double> qual_sd;
    std::vector<double> mismatch_qual_mean;
    std::vector<double> mismatch_qual_sd;

    SeqValues(const SeqOptions& opts)
        : read_length(opts.read_length), mismatch_probs(opts.read_length),
          prob_insert(opts.prob_insert), prob_delete(opts.prob_delete),
          qual_mean(), qual_sd(),
          mismatch_qual_mean(), mismatch_qual_sd() {

        // Compute probability at raise point.
        double y_r = 2 * opts.prob_mismatch -
            opts.position_raise * opts.prob_mismatch_begin -
            opts.prob_mismatch_end +
            opts.prob_mismatch_end * opts.position_raise;
        /*
         Compute mismatch probability at each base.
         Use piecewise linear function for mismatch probability simulation.
         */
        for (uint32 i = 0; i < read_length; i++) {
            double x = static_cast<double>(i) / (opts.read_length - 1);
            if (x < opts.position_raise) {
                double b = opts.prob_mismatch_begin;
                double m = (y_r - opts.prob_mismatch_begin) /
                    opts.position_raise;
                mismatch_probs[i] = m * x + b;
            } else {
                double b = y_r;
                double m = (opts.prob_mismatch_end - y_r) /
                    (1 - opts.position_raise);
                x -= opts.position_raise;
                mismatch_probs[i] = m * x + b;
            }
        }
        if (opts.prob_mismatch_scale != 1.0) {
            for (uint32 i = 0; i < read_length; ++i) {
                mismatch_probs[i] *= opts.prob_mismatch_scale;
            }
        }
        // Compute match/mismatch means and standard deviations.
        compute_means_sds(mismatch_qual_mean,
                          opts.mean_mismatch_qual_begin,
                          opts.mean_mismatch_qual_end);
        compute_means_sds(mismatch_qual_sd,
                          opts.sd_mismatch_qual_begin,
                          opts.sd_mismatch_qual_end);
        compute_means_sds(qual_mean, opts.mean_qual_begin, opts.mean_qual_end);
        compute_means_sds(qual_sd, opts.sd_qual_begin, opts.sd_qual_end);

    }


private:

    void compute_means_sds(std::vector<double>& to_change,
                           const double& x_begin,
                           const double& x_end) {
        to_change.reserve(read_length);
        for (uint32 i = 0; i < read_length; ++i) {
            double b = x_begin;
            double x = static_cast<double>(i) / static_cast<double>(read_length - 1);
            double m = (x_end - x_begin);
            to_change.push_back(m * x + b);
        }
    }

};




/*
 Way to store info for the sequence identifier line in FASTQ file.
 It's organized as such:
 @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
    <read>:<is filtered>:<control number>:<sample number>
 (where the newline and tab are, just have one space)
 (from https://help.basespace.illumina.com/articles/descriptive/fastq-files/)
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

    uint32 sample(pcg32& eng) const {
        uint32 rnd = sampler.sample(eng);
        uint32 len_ = region_len * rnd;
        if (region_len > 1) len_ += static_cast<uint32>(runif_01(eng) * region_len);
        return len_;
    }

};



// for Illumina read lengths: return constant value unless fragment length < constant
struct IlluminaReads {

    uint32 read_length;         // read length (constant)

    IlluminaReads() {};
    IlluminaReads(const uint32& read_length_) : read_length(read_length_) {};
    // Copy constructor
    IlluminaReads(const IlluminaReads& other)
        : read_length(other.read_length) {};
    // Assignment operator
    IlluminaReads& operator=(const IlluminaReads& other) {
        read_length = other.read_length;
        return *this;
    }

    inline uint32 sample(const uint32& frag_len) const {
        if (frag_len < read_length) return frag_len;
        return read_length;
    }

};


// for long-read sequencing (PacBio or Nanopore) read lengths: simply return fragment length
struct LongReadReads {

    LongReadReads() {};
    LongReadReads(const uint32& read_length_) {};
    // Copy constructor
    LongReadReads(const LongReadReads& other) {};
    // Assignment operator
    LongReadReads& operator=(const LongReadReads& other) {
        return *this;
    }

    inline uint32 sample(const uint32& frag_len) const {
        return frag_len;
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
        : region_len(region_len_), means(means_), sds(sds_), distr(0, 1) {
        if (means_.size() != sds_.size()) stop("means and sds must be the same length.");
    }
    QualitySampler(const List& qual_info)
        : region_len(), means(), sds(), distr(0, 1) {
        const std::vector<double>& means_(qual_info["means"]);
        const std::vector<double>& sds_(qual_info["sds"]);
        const uint32& region_len_(qual_info["region_len"]);
        if (means_.size() != sds_.size()) stop("means and sds must be the same length.");
        QualitySampler other(means_, sds_, region_len_);
        *this = other;
    }
    QualitySampler& operator=(const QualitySampler& other) {
        region_len = other.region_len;
        means = other.means;
        sds = other.sds;
        return *this;
    }



    char sample(uint32 pos, pcg32& eng) const {
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
        // Limit quality to [0, 40]
        if (rnd < 0) rnd = 0;
        if (rnd > 40) rnd = 40;
        // Convert to int, then char
        char q = static_cast<char>('!' + static_cast<int>(rnd));
        return q;
    }

private:

    mutable std::normal_distribution<double> distr;

};






/*
 Sample for sequencing error when rates vary by position on read.
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
    SequencingErrorSampler(const std::vector<double>& match_probs,
                           const std::vector<double>& mis_probs,
                           const std::vector<double>& ins_probs,
                           const std::vector<double>& del_probs,
                           const uint32& region_len_)
        : region_len(region_len_), probs() {
        if (match_probs.size() != mis_probs.size() ||
            match_probs.size() != ins_probs.size() ||
            match_probs.size() != del_probs.size()) {
            stop("all input probs must be the same length.");
        }

        probs.reserve(match_probs.size());
        std::vector<double> new_probs(3);

        for (uint32 i = 0; i < mis_probs.size(); i++) {
            double mat_i = match_probs[i];
            double mis_i = mis_probs[i];
            double ins_i = ins_probs[i];
            double del_i = del_probs[i];
            /*
             Making sure they sum to 1 (deletion don't need to be corrected bc they're
             not included in probs, they're implicitly known once each set of
             probs sum to 1):
             */
            double sum_i = mat_i + mis_i + ins_i + del_i;
            mat_i /= sum_i;
            mis_i /= sum_i;
            ins_i /= sum_i;
            new_probs[0] = mat_i;
            new_probs[1] = mat_i + mis_i;
            new_probs[2] = mat_i + mis_i + ins_i;
            probs.push_back(new_probs);
        }

    }

    SequencingErrorSampler(const List& seq_error_info) {
        const std::vector<double>& match_probs(seq_error_info["match_probs"]);
        const std::vector<double>& mis_probs(seq_error_info["mis_probs"]);
        const std::vector<double>& ins_probs(seq_error_info["ins_probs"]);
        const std::vector<double>& del_probs(seq_error_info["del_probs"]);
        const uint32& region_len_(seq_error_info["region_len"]);
        SequencingErrorSampler other(match_probs, mis_probs, ins_probs, del_probs,
                                     region_len_);
        *this = other;
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
    inline uint8 sample(uint32 pos, pcg32& eng) const {

        // Adjust pos to point to inside `probs` vector
        pos /= region_len;
        if (pos >= probs.size()) {
            Rcout << "`pos` is past SequencingErrorSampler bounds; using last value";
            Rcout << std::endl;
            pos = probs.size() - 1;
        }

        uint8 out;
        double u = runif_01(eng);
        if (u < probs[pos][1]) { // match
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

    WGS_t() {};
    WGS_t(const T& seq_object,
          const std::vector<double>& frag_len_probs,
          const uint32& frag_len_region_len,
          const List& seq_error_info,
          const std::vector<double>& ins_length_probs,
          const std::vector<double>& del_length_probs,
          const List& qual_info,
          const List& mis_qual_info)
        : seqs(),
          frag_lengths(frag_len_probs, frag_len_region_len),
          errors(seq_error_info),
          ins_lengths(ins_length_probs),
          del_lengths(del_length_probs),
          quals(qual_info),
          mis_quals(mis_qual_info),
          seq_lengths(seq_object.seq_sizes()),
          sequences(&seq_object),
          mm(256, "") {

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

    // Sample a starting location for a fragment
    inline uint32 frag_loc_sample(const uint32& frag_len,
                                  const uint32& seq_len,
                                  pcg32& eng) {
        if (frag_len >= seq_len) return 0U;
        double u = runif_01(eng);
        uint32 pos = static_cast<uint32>(u * (seq_len - frag_len + 1));
        return pos;
    }

    // Adjust one character for a mismatch
    void mismatch(char& input, pcg32& eng) {
        const std::string& str(mm[static_cast<uint32>(input)]);
        uint32 ind = static_cast<uint32>(runif_01(eng) * str.size());
        input = str[ind];
        return;
    }

    // Adjust string and qualities for an insertion
    void insertion(const sint32& size, const uint32& read_pos,
                   std::string& read, std::string& qual,
                   pcg32& eng) {
        uint32 size_ = static_cast<uint32>(size);
        std::string nts;
        nts.reserve(size_);
        for (uint32 j = 0; j < size_; j++) {
            uint32 k = runif_01(eng) * 4;
            nts += sequencer::bases[k];
            qual += this->mis_quals.sample(read_pos, eng);
        }
        read.insert(read_pos + 1, nts);
        return;
    }

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
                     pcg32& eng) {

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
                  const List& seq_error_info,
                  const std::vector<double>& ins_length_probs,
                  const std::vector<double>& del_length_probs,
                  const List& qual_info,
                  const List& mis_qual_info,
                  const uint32& read_length_,
                  const bool& paired)
        : WGS_t<T>(
                seq_object, frag_len_probs,
                frag_len_region_len, seq_error_info, ins_length_probs, del_length_probs,
                qual_info, mis_qual_info), read_length(read_length_), n_reads(1) {
                    if (paired) n_reads = 2;
                };
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


    // Sample one set of read strings (4 lines: ID, sequence, "+", quality)
    std::vector<std::string> sample(pcg32& eng, SequenceIdentifierInfo& ID_info) {

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
                   pcg32& eng) {

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
                           pcg32& eng) {
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


typedef IlluminaWGS_t<VarGenome> VariantIlluminaWGS;
typedef IlluminaWGS_t<RefGenome> ReferenceIlluminaWGS;




















#endif
