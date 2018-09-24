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
 Sample for lengths that can be quite long.
 This is used to avoid having to store a probability for every single potential length,
 if the lengths can extend to something like 1 Mbp.
 This class is used to sample fragment lengths and Nanopore/PacBio read lengths.

 It does weighted sampling using a vector of probabilities for each bin of lengths.
 Once it chooses a bin, it then uniformly samples within that bin if bin length > 1.
 */
class LongLenSampler {

public:

    uint32 region_len;                  // length of each region
    TableSampler sampler;               // sampler that chooses a region

    LongLenSampler() {};
    LongLenSampler(const std::vector<double>& probs,
                    const uint32& region_len_)
        : region_len(region_len_), sampler(probs) {}
    // Copy constructor
    LongLenSampler(const LongLenSampler& other)
        : region_len(other.region_len), sampler(other.sampler) {};
    // Assignment operator
    LongLenSampler& operator=(const LongLenSampler& other) {
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

typedef LongLenSampler FragLenSampler;  // for fragment lengths
typedef LongLenSampler LongReadSampler; // for PacBio or Nanopore read lengths




class IlluminaReadSampler {

public:

    uint32 read_length;         // read length (constant)

    IlluminaReadSampler() {};
    IlluminaReadSampler(const uint32& read_length_) : read_length(read_length_) {};
    // Copy constructor
    IlluminaReadSampler(const IlluminaReadSampler& other)
        : read_length(other.read_length) {};
    // Assignment operator
    IlluminaReadSampler& operator=(const IlluminaReadSampler& other) {
        read_length = other.read_length;
        return *this;
    }

    inline uint32 sample() const {
        return read_length;
    }
    inline uint32 sample(pcg32& eng) const {
        return read_length;
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

    QualitySampler() : region_len(), means(), sds(), distr(1.0, 0) {};
    QualitySampler(const std::vector<double>& means_,
                   const std::vector<double>& sds_,
                   const uint32& region_len_)
        : region_len(region_len_), means(means_), sds(sds_), distr(1.0, 0) {
        if (means_.size() != sds_.size()) stop("means and sds must be the same length.");
    }
    QualitySampler(const List& qual_info)
        : region_len(), means(), sds(), distr(1.0, 0) {
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
     Returns 0 for match, 1 for mismatch, 2 for insertion, 3 for deletion
     */
    inline uint32 sample(uint32 pos, pcg32& eng) const {

        uint32 out;

        // Adjust pos to point to inside `probs` vector
        pos /= region_len;
        if (pos >= probs.size()) {
            Rcout << "`pos` is past SequencingErrorSampler bounds; using last value";
            Rcout << std::endl;
            pos = probs.size() - 1;
        }
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
 Class to do everything for whole-genome sequencing.

 `T` should be `IlluminaReadSampler` or `LongReadSampler`.

 `U` should be `VarGenome` or `RefGenome`
 */
template <typename T, typename U>
class WholeGenomeSequencerType {
public:

    // Samplers:                       Returns:
    T read_lengths;                 // read lengths
    TableSampler seqs;              // index for which genome-sequence to sequence
    FragLenSampler frag_lengths;    // fragment lengths
    SequencingErrorSampler errors;  // sequencing errors
    TableSampler ins_sizes;         // insertion sizes
    TableSampler del_sizes;         // deletion sizes
    QualitySampler quals;           // qualities for matches
    QualitySampler mis_quals;       // qualities for mismatches

    // Info:
    std::vector<uint32> seq_sizes;  // genome-sequence sizes
    const U* sequences;             // pointer to `const U`//
    std::vector<std::string> mm;    // mismatch vector

    WholeGenomeSequencerType() {};
    WholeGenomeSequencerType(const U& seq_object,
                             const std::vector<double>& read_len_probs,
                             const uint32& read_len_region_len,
                             const std::vector<double>& frag_len_probs,
                             const uint32& frag_len_region_len,
                             const List& seq_error_info,
                             const std::vector<double>& ins_size_probs,
                             const std::vector<double>& del_size_probs,
                             const List& qual_info,
                             const List& mis_qual_info)
        : read_lengths(read_len_probs, read_len_region_len),
          seqs(),
          frag_lengths(frag_len_probs, frag_len_region_len),
          errors(seq_error_info),
          ins_sizes(ins_size_probs),
          del_sizes(del_size_probs),
          quals(qual_info),
          mis_quals(mis_qual_info),
          seq_sizes(seq_object.seq_sizes()),
          sequences(&seq_object),
          mm(256, "") {
        std::vector<double> probs = static_cast<std::vector<double>>(seq_sizes);
        seqs = TableSampler(probs);

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

    void mismatch(char& input, pcg32& eng) {
        const std::string& str(mm[static_cast<uint32>(input)]);
        uint32 ind = static_cast<uint32>(runif_01(eng) * str.size());
        input = str[ind];
        return;
    }




};

template <typename T>
class VarWholeGenomeSequencer: public WholeGenomeSequencerType<T, VarGenome> {
public:

    // Constructors:
    VarWholeGenomeSequencer() : WholeGenomeSequencerType<T, VarGenome>() {}
    VarWholeGenomeSequencer(const MutationRates& mr_)
        : WholeGenomeSequencerType<T, VarGenome>(mr_, 0) {};
    // Copy constructor
    VarWholeGenomeSequencer(const VarWholeGenomeSequencer<T>& other)
        : WholeGenomeSequencerType<T, VarGenome>(other) {};
    // Assignment operator
    VarWholeGenomeSequencer<T>& operator=(const VarWholeGenomeSequencer<T>& other) {
        WholeGenomeSequencerType<T, VarGenome>::operator=(other);
        return *this;
    }


    void sample(std::string& read, std::string& qual, pcg32& eng) {

        uint32 read_len = this->read_lengths.sample(eng);
        uint32 seq_ind = this->seqs.sample(eng);
        uint32 frag_len = this->frag_lengths.sample(eng);
        const VarSequence& var_seq((*(this->sequences))[seq_ind]);
        uint32 seq_len = var_seq.size();
        // Clear read and quality strings, in case they were previously altered:
        read.clear();
        qual.clear();
        // Now fill `read` from `var_seq`:
        uint32 mut_i = 0;
        uint32 start = frag_loc_sample(frag_len, this->seq_sizes[seq_ind], eng);
        var_seq.set_seq_chunk(read, start, read_len, mut_i);
        // Adjust in case it's at the end of the sequence:
        read_len = read.size();
        // Reserve memory for quality:
        qual.reserve(read_len);
        // Position on the `var_seq` sequence of the last nucleotide in `read`:
        uint32 seq_pos = start + read_len - 1;
        // For current position on `read`:
        uint32 i = 0;
        while (i < read_len && seq_pos < seq_len) {
            uint32 err_code = this->errors.sample(i, eng);
            if (err_code == 0) { // -------------------- match
                qual += this->quals.sample(i, eng);
                i++;
                seq_pos++;
            } else if (err_code == 1) { // ------------- mismatch
                this->mismatch(read[i], eng);
                qual += this->mis_quals.sample(i, eng);
                i++;
                seq_pos++;
            } else if (err_code == 2) { // ------------- insertion
                uint32 size_ = this->ins_sizes(eng) + 1;
                std::string nts;
                nts.reserve(size_);
                for (uint32 j = 0; j < size_; j++) {
                    uint32 k = runif_01(eng) * 4;
                    nts += sequencer::bases[k];
                    qual += this->mis_quals.sample(i, eng);
                }
                read.resize(read_len - size_);  // <-- so `read` remains the same size
                read.insert(i + 1, nts);
                i += (size_ + 1);
                seq_pos++;
            } else if (err_code == 3) { // ------------- deletion
                uint32 size_ = this->del_sizes(eng) + 1;
                if (size_ > (seq_len - i)) size_ = seq_len - i;
                read.erase(i, size_);
                // char get_nt(const uint32& new_pos)
                // So `read` remains the same size:
                for (uint32 j = 0; j < size_; j++) {
                    if (seq_pos >= seq_len) break;
                    read += var_seq.get_nt(seq_pos);
                    seq_pos++;
                }
            }
        }

    }
};


















#endif
