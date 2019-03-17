

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
#include "pacbio.hpp"  // PacBio* types




uint32 PacBioReadLenSampler::sample(pcg64& eng) {
    uint32 len_;
    if (use_distr) {
        double rnd = distr(eng) + loc;
        uint32 iters = 0; // to make sure it doesn't run many times
        // Rejection sampling to keep it above minimum and not NaN:
        while (rnd < min_read_len && iters < 10) {
            rnd = distr(eng) + loc;
            iters++;
        }
        // Give up if it keeps returning craziness:
        if (rnd < min_read_len) rnd = min_read_len;
        len_ = static_cast<uint32>(rnd);
    } else {
        uint64 ind = sampler.sample(eng);
        len_ = read_lens[ind];
    }
    return len_;
}



double PacBioQualityError::calc_min_exp() {

    double min_exp_ = 1;

    double total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
        std::pow(prob_subst, min_exp_);
    double left, right;

    if (total < prob_thresh) {
        while (total < prob_thresh) {
            min_exp_ /= 2;
            total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
                std::pow(prob_subst, min_exp_);
        }
        left = min_exp_;
        right = min_exp_ * 2;
    } else {
        while (total > prob_thresh) {
            min_exp_ *= 2;
            total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
                std::pow(prob_subst, min_exp_);
        }
        left = min_exp_ / 2;
        right = min_exp_;
    }
    for (uint32 i = 0; i < 15; i++) {
        double m = (left + right) / 2;
        total = std::pow(prob_ins, m) + std::pow(prob_del, m) + std::pow(prob_subst, m);
        if (total == prob_thresh){
            min_exp_ = m;
            break;
        } else if (total > prob_thresh) {
            left = m;
            min_exp_ = (m + right) / 2;
        } else {
            right = m;
            min_exp_ = (left + m) / 2;
        }
    }

    return min_exp_;
}




void PacBioQualityError::update_probs(pcg64& eng,
                                      const double& passes_left,
                                      const double& passes_right) {


    // Thresholds for truncated normal distribution
    double left_thresh = (min_exp - (std::sqrt(passes_left + sqrt_params[0]) -
                          sqrt_params[1]))  / sigmoid(passes_left);
    double right_thresh = (min_exp - (std::sqrt(passes_right + sqrt_params[0]) -
                           sqrt_params[1])) / sigmoid(passes_right);

    // Sample for noise:
    double incr_quals_l = trunc_norm(left_thresh, eng);
    double incr_quals_r = trunc_norm(right_thresh, eng);

    // Exponents of quality increase
    double exponent_left = incr_quals_l * sigmoid(passes_left) +
        std::sqrt(passes_left + sqrt_params[0]) - sqrt_params[1];
    double exponent_right = incr_quals_r * sigmoid(passes_right) +
        std::sqrt(passes_right + sqrt_params[0]) - sqrt_params[1];
    // SimLoRD prevents from going below 0.6 by chance:
    if (exponent_left < 0.6) exponent_left = 0.6;
    if (exponent_right < 0.6) exponent_right = 0.6;

    // Now calculate new cumulative error probabilities:
    cum_probs_left[0] = std::pow(prob_ins, exponent_left);
    cum_probs_left[1] = std::pow(prob_del, exponent_left) + cum_probs_left[0];
    cum_probs_left[2] = std::pow(prob_subst, exponent_left) + cum_probs_left[1];

    cum_probs_right[0] = std::pow(prob_ins, exponent_right);
    cum_probs_right[1] = std::pow(prob_del, exponent_right) + cum_probs_right[0];
    cum_probs_right[2] = std::pow(prob_subst, exponent_right) + cum_probs_right[1];

    return;

}


template <typename T>
void PacBioOneGenome<T>::one_read(std::vector<std::string>& fastq_chunks,
                                  pcg64& eng,
                                  SequenceIdentifierInfo& ID_info) {

    std::string& fastq_chunk(fastq_chunks[0]);

    /*
    Sample read info, and set the sequence space(s) required for these read(s).
    */
    // Sample sequence:
    seq_ind = seq_sampler.sample(eng);
    uint32 seq_len = (*sequences)[seq_ind].size();

    // Sample read length:
    read_length = len_sampler.sample(eng);
    if (read_length >= seq_len) read_length = seq_len;

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      seq_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/variant sequence needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more sequence bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_seq_space = read_length + deletions.size() - insertions.size();

    // Sample read starting position:
    if (read_seq_space < seq_len) {
        double u = runif_01(eng);
        read_start = static_cast<uint32>(u * (seq_len - read_seq_space + 1));
    } else if (read_seq_space == seq_len) {
        read_start = 0;
    } else {
        stop("read_seq_space should never exceed the sequence length.");
    }

    // Fill the reads and qualities
    append_chunk(fastq_chunk, eng, ID_info);

    return;
}



template <typename T>
void PacBioOneGenome<T>::re_read(std::vector<std::string>& fastq_chunks,
                                 pcg64& eng,
                                 SequenceIdentifierInfo& ID_info) {

    std::string& fastq_chunk(fastq_chunks[0]);

    /*
     Use the same read info as before.
    */
    uint32 seq_len = (*sequences)[seq_ind].size();

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      seq_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/variant sequence needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more sequence bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_seq_space = read_length + deletions.size() - insertions.size();

    /*
     In the very rare situation where a duplication occurs, then enough deletions
     happen where the required sequence space exceeds what's available, I'm going
     to remove deletions until we have enough room.
     */
    while ((read_seq_space + read_start) > seq_len) {
        if (deletions.empty()) break;
        deletions.pop_back();
        read_seq_space--;
    }
    // If that still doesn't work, I give up on the duplicate.
    if ((read_seq_space + read_start) > seq_len) return;

    // Fill the reads and qualities
    append_chunk(fastq_chunk, eng, ID_info);

    return;
}




template <typename T>
void PacBioOneGenome<T>::append_chunk(std::string& fastq_chunk,
                                      pcg64& eng,
                                      SequenceIdentifierInfo& ID_info) {

    // Make sure it has enough memory reserved:
    fastq_chunk.reserve(fastq_chunk.size() + read_length * 3 + 10);

    fastq_chunk += ID_info.get_line() + '\n';

    // Boolean for whether we take the reverse side:
    bool reverse = runif_01(eng) < 0.5;

    // Fill in read:
    (*sequences)[seq_ind].fill_read(read, 0, read_start, read_seq_space);

    // Reverse complement if necessary:
    if (reverse) rev_comp(read, read_seq_space);

    /*
     Adding read with errors:
     */
    uint32 read_pos = 0;
    uint32 current_length = 0;
    uint32 rndi;
    while (current_length < read_length) {
        if (!insertions.empty() && read_pos == insertions.front()) {
            rndi = runif_aabb(eng, static_cast<uint32>(0UL), static_cast<uint32>(3UL));
            fastq_chunk += read[read_pos];
            fastq_chunk += alias_sampler::bases[rndi];
            insertions.pop_front();
            current_length += 2;
        } else if (!deletions.empty() && read_pos == deletions.front()) {
            deletions.pop_front();
        } else if (!substitutions.empty() && read_pos == substitutions.front()) {
            rndi = runif_aabb(eng, static_cast<uint32>(0UL), static_cast<uint32>(2UL));
            fastq_chunk += mm_nucleos[nt_map[read[read_pos]]][rndi];
            substitutions.pop_front();
            current_length++;
        } else {
            fastq_chunk += read[read_pos];
            current_length++;
        }
        read_pos++;
    }

    fastq_chunk += "\n+\n";

    // Adding qualities:
    for (uint i = 0; i < split_pos; i++) fastq_chunk += qual_left;
    for (uint i = split_pos; i < read_length; i++) fastq_chunk += qual_right;
    fastq_chunk += '\n';

    return;
}






/*
 ========================================================================================
 ========================================================================================

 Writing reads

 ========================================================================================
 ========================================================================================
 */




//' PacBio sequence for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void pacbio_ref_cpp(SEXP ref_genome_ptr,
                      const std::string& out_prefix,
                      const bool& compress,
                      const uint32& n_reads,
                      const uint32& n_cores,
                      const bool& show_progress,
                      const uint32& read_chunk_size,
                      const double& prob_dup,
                      const double& scale,
                      const double& sigma,
                      const double& loc,
                      const double& min_read_len,
                      const std::vector<double>& read_probs,
                      const std::vector<uint32>& read_lens,
                      const uint32& max_passes,
                      const std::vector<double>& chi2_params_n,
                      const std::vector<double>& chi2_params_s,
                      const std::vector<double>& sqrt_params,
                      const std::vector<double>& norm_params,
                      const double& prob_thresh,
                      const double& prob_ins,
                      const double& prob_del,
                      const double& prob_subst,
                      const std::string& instrument,
                      const uint32& run_number,
                      const std::string& flowcell_ID,
                      const uint32& lane,
                      const uint32& tile,
                      const uint32& x_pos,
                      const uint32& y_pos,
                      const uint32& read,
                      const std::string& is_filtered,
                      const uint32& control_number,
                      const uint32& sample_number) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    PacBioReference read_filler_base;

    if (read_probs.size() == 0) {
        read_filler_base =
            PacBioReference(*ref_genome,
                            scale, sigma, loc, min_read_len,
                            max_passes, chi2_params_n, chi2_params_s,
                            sqrt_params, norm_params, prob_thresh,
                            prob_ins, prob_del, prob_subst);
    } else {
        read_filler_base =
            PacBioReference(*ref_genome,
                            read_probs, read_lens,
                            max_passes, chi2_params_n, chi2_params_s,
                            sqrt_params, norm_params, prob_thresh,
                            prob_ins, prob_del, prob_subst);
    }

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);

    if (compress) {
        write_reads_cpp_<PacBioReference, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads,
                prob_dup, read_chunk_size, 1U, n_cores, show_progress);
    } else {
        write_reads_cpp_<PacBioReference, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads,
                prob_dup, read_chunk_size, 1U, n_cores, show_progress);
    }


    return;
}




//' PacBio sequence for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void pacbio_var_cpp(SEXP var_set_ptr,
                    const std::string& out_prefix,
                    const bool& compress,
                    const uint32& n_reads,
                    const uint32& n_cores,
                    const bool& show_progress,
                    const uint32& read_chunk_size,
                    const std::vector<double>& variant_probs,
                    const double& prob_dup,
                    const double& scale,
                    const double& sigma,
                    const double& loc,
                    const double& min_read_len,
                    const std::vector<double>& read_probs,
                    const std::vector<uint32>& read_lens,
                    const uint32& max_passes,
                    const std::vector<double>& chi2_params_n,
                    const std::vector<double>& chi2_params_s,
                    const std::vector<double>& sqrt_params,
                    const std::vector<double>& norm_params,
                    const double& prob_thresh,
                    const double& prob_ins,
                    const double& prob_del,
                    const double& prob_subst,
                    const std::string& instrument,
                    const uint32& run_number,
                    const std::string& flowcell_ID,
                    const uint32& lane,
                    const uint32& tile,
                    const uint32& x_pos,
                    const uint32& y_pos,
                    const uint32& read,
                    const std::string& is_filtered,
                    const uint32& control_number,
                    const uint32& sample_number) {

    XPtr<VarSet> var_set(var_set_ptr);
    PacBioVariants read_filler_base;

    if (read_probs.size() == 0) {
        read_filler_base =
            PacBioVariants(*var_set, variant_probs,
                           scale, sigma, loc, min_read_len,
                           max_passes, chi2_params_n, chi2_params_s,
                           sqrt_params, norm_params, prob_thresh,
                           prob_ins, prob_del, prob_subst);
    } else {
        read_filler_base =
            PacBioVariants(*var_set, variant_probs,
                           read_probs, read_lens,
                           max_passes, chi2_params_n, chi2_params_s,
                           sqrt_params, norm_params, prob_thresh,
                           prob_ins, prob_del, prob_subst);
    }

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);

    if (compress) {
        write_reads_cpp_<PacBioVariants, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads,
                prob_dup, read_chunk_size, 1U, n_cores, show_progress);
    } else {
        write_reads_cpp_<PacBioVariants, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads,
                prob_dup, read_chunk_size, 1U, n_cores, show_progress);
    }

    return;
}
