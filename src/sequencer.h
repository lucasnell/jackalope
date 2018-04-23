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

 Edited for use in gemino by Lucas Nell, 2017

*/

#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng
#include <unordered_map> // unordered_map
#include <string>  // string class

#include "gemino_types.h"  // uint

using namespace Rcpp;


namespace sequencer
{
    std::string bases = "ACGTN";
    // For choosing a base for a mismatch:
    std::unordered_map<char,std::string> pick_mismatch = {
        {'A', "CGTN"}, {'C', "AGTN"}, {'G', "ACTN"}, {'T', "ACGN"}, {'N', "ACGT"}
    };
}


struct seq_options
{
    // Constructor
    seq_options() {

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

};


struct seq_values
{
    seq_values(seq_options opts) {
        read_length = static_cast<uint>(opts.read_length);
        prob_insert = opts.prob_insert;
        prob_delete = opts.prob_delete;
        mismatch_probs = std::vector<double>(read_length);
        qual_mean = std::vector<double>(read_length);
        qual_sd = std::vector<double>(read_length);
        mismatch_qual_mean = std::vector<double>(read_length);
        mismatch_qual_sd = std::vector<double>(read_length);

        // Compute probability at raise point.
        double y_r = 2 * opts.prob_mismatch -
            opts.position_raise * opts.prob_mismatch_begin -
            opts.prob_mismatch_end +
            opts.prob_mismatch_end * opts.position_raise;
        // Compute mismatch probability at each base.
        // Use piecewise linear function for mismatch probability simulation.
        for (uint i = 0; i < read_length; i++) {
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
            for (uint i = 0; i < read_length; ++i) {
                mismatch_probs[i] *= opts.prob_mismatch_scale;
            }
        }
        // Compute match/mismatch means and standard deviations.
        for (uint i = 0; i < read_length; ++i) {
            double b = opts.mean_mismatch_qual_begin;
            double x = static_cast<double>(i) / (opts.read_length - 1);
            double m = (opts.mean_mismatch_qual_end -
                        opts.mean_mismatch_qual_begin);
            mismatch_qual_mean[i] = m * x + b;
        }
        for (uint i = 0; i < read_length; ++i) {
            double b = opts.sd_mismatch_qual_begin;
            double x = static_cast<double>(i) / (opts.read_length - 1);
            double m = (opts.sd_mismatch_qual_end -
                        opts.sd_mismatch_qual_begin);
            mismatch_qual_sd[i] = m * x + b;
        }
        for (uint i = 0; i < read_length; ++i) {
            double b = opts.mean_qual_begin;
            double x = static_cast<double>(i) / (opts.read_length - 1);
            double m = (opts.mean_qual_end -
                        opts.mean_qual_begin);
            qual_mean[i] = m * x + b;
        }
        for (uint i = 0; i < read_length; ++i) {
            double b = opts.sd_qual_begin;
            double x = static_cast<double>(i) / (opts.read_length - 1);
            double m = (opts.sd_qual_end -
                        opts.sd_qual_begin);
            qual_sd[i] = m * x + b;
        }
    }

    uint read_length;
    std::vector<double> mismatch_probs;
    double prob_insert;
    double prob_delete;
    std::vector<double> qual_mean;
    std::vector<double> qual_sd;
    std::vector<double> mismatch_qual_mean;
    std::vector<double> mismatch_qual_sd;

};



#endif
