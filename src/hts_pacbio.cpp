
#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng
#include <string>  // string class
#include <random>  // distributions



#include "jackalope_types.h"  // uint64
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // runif_01
#include "alias_sampler.h"  // AliasSampler
#include "util.h"  // clear_memory
#include "str_manip.h"  // rev_comp
#include "hts.h"  // generic sequencer classes
#include "hts_pacbio.h"  // PacBio* types
#include "io.h"  // File* types




uint64 PacBioReadLenSampler::sample(pcg64& eng) {
    uint64 len_;
    if (use_distr) {
        double rnd = distr(eng) + loc;
        uint64 iters = 0; // to make sure it doesn't run many times
        // Rejection sampling to keep it above minimum and not NaN:
        while (rnd < min_read_len && iters < 10) {
            rnd = distr(eng) + loc;
            iters++;
        }
        // Give up if it keeps returning craziness:
        if (rnd < min_read_len) rnd = min_read_len;
        len_ = static_cast<uint64>(rnd);
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
    for (uint64 i = 0; i < 15; i++) {
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
template <typename U>
void PacBioOneGenome<T>::one_read(std::vector<U>& fastq_pools,
                                  bool& finished,
                                  pcg64& eng) {

    U& fastq_pool(fastq_pools[0]);

    /*
    Sample read info, and set the chromosome space(s) required for these read(s).
    */
    // Get chromosome:
    chrom_ind = 0;
    while (chrom_ind < chrom_reads.size() && chrom_reads[chrom_ind] == 0) chrom_ind++;
    if (chrom_ind == chrom_reads.size()) {
        finished = true;
        return;
    }

    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample read length:
    read_length = len_sampler.sample(eng);
    if (read_length >= chrom_len) read_length = chrom_len;

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      chrom_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/haplotype chromosome needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more chromosome bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_chrom_space = read_length + deletions.size() - insertions.size();

    // Sample read starting position:
    if (read_chrom_space < chrom_len) {
        double u = runif_01(eng);
        read_start = static_cast<uint64>(u * (chrom_len - read_chrom_space + 1));
    } else if (read_chrom_space == chrom_len) {
        read_start = 0;
    } else {
        stop("read_chrom_space should never exceed the chromosome length.");
    }

    // Fill the reads and qualities
    append_pool<U>(fastq_pool, eng);

    return;
}

// Overloaded for when we input a haplotype chromosome stored as string
template <typename T>
template <typename U>
void PacBioOneGenome<T>::one_read(const std::string& chrom,
                                  const uint64& chrom_i,
                                  std::vector<U>& fastq_pools,
                                  pcg64& eng) {

    U& fastq_pool(fastq_pools[0]);

    chrom_ind = chrom_i;

    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample read length:
    read_length = len_sampler.sample(eng);
    if (read_length >= chrom_len) read_length = chrom_len;

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      chrom_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/haplotype chromosome needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more chromosome bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_chrom_space = read_length + deletions.size() - insertions.size();

    // Sample read starting position:
    if (read_chrom_space < chrom_len) {
        double u = runif_01(eng);
        read_start = static_cast<uint64>(u * (chrom_len - read_chrom_space + 1));
    } else if (read_chrom_space == chrom_len) {
        read_start = 0;
    } else {
        stop("read_chrom_space should never exceed the chromosome length.");
    }

    // Fill the reads and qualities
    append_pool<U>(chrom, fastq_pool, eng);

    return;
}






template <typename T>
template <typename U>
void PacBioOneGenome<T>::re_read(std::vector<U>& fastq_pools,
                                 bool& finished,
                                 pcg64& eng) {

    U& fastq_pool(fastq_pools[0]);

    /*
     Use the same read info as before.
    */
    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      chrom_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/haplotype chromosome needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more chromosome bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_chrom_space = read_length + deletions.size() - insertions.size();

    /*
     In the very rare situation where a duplication occurs, then enough deletions
     happen where the required chromosome space exceeds what's available, I'm going
     to remove deletions until we have enough room.
     */
    while ((read_chrom_space + read_start) > chrom_len) {
        if (deletions.empty()) break;
        deletions.pop_back();
        read_chrom_space--;
    }
    // If that still doesn't work, I give up on the duplicate.
    if ((read_chrom_space + read_start) > chrom_len) return;

    // Fill the reads and qualities
    append_pool<U>(fastq_pool, eng);

    return;
}




template <typename T>
template <typename U>
void PacBioOneGenome<T>::re_read(const std::string& chrom,
                                 const uint64& chrom_i,
                                 std::vector<U>& fastq_pools,
                                 pcg64& eng) {

    U& fastq_pool(fastq_pools[0]);

    chrom_ind = chrom_i;

    /*
     Use the same read info as before.
     */
    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample for # passes over read:
    pass_sampler.sample(split_pos, passes_left, passes_right, eng, read_length);

    // Sample for errors and qualities:
    qe_sampler.sample(eng, qual_left, qual_right, insertions, deletions, substitutions,
                      chrom_len, read_length, split_pos, passes_left, passes_right);

    /*
     The amount of space on the reference/haplotype chromosome needed to create this read.
     I'm adding deletions because more deletions mean that I need
     more chromosome bases to achieve the same read length.
     Insertions means I need fewer.
     */
    read_chrom_space = read_length + deletions.size() - insertions.size();

    /*
     In the very rare situation where a duplication occurs, then enough deletions
     happen where the required chromosome space exceeds what's available, I'm going
     to remove deletions until we have enough room.
     */
    while ((read_chrom_space + read_start) > chrom_len) {
        if (deletions.empty()) break;
        deletions.pop_back();
        read_chrom_space--;
    }
    // If that still doesn't work, I give up on the duplicate.
    if ((read_chrom_space + read_start) > chrom_len) return;

    // Fill the reads and qualities
    append_pool<U>(chrom, fastq_pool, eng);

    return;
}





template <typename T>
template <typename U>
void PacBioOneGenome<T>::append_pool(U& fastq_pool, pcg64& eng) {

    // Make sure it has enough memory reserved:
    fastq_pool.reserve(fastq_pool.size() + read_length * 3 + 10);

    // Boolean for whether we take the reverse side:
    bool reverse = runif_01(eng) < 0.5;

    // ID line:
    fastq_pool.push_back('@');
    for (const char& c : name) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    for (const char& c : (*chromosomes)[chrom_ind].name) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    for (const char& c : std::to_string(read_start)) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    if (reverse) {
        fastq_pool.push_back('R');
    } else fastq_pool.push_back('F');
    fastq_pool.push_back('\n');

    // Fill in read:
    (*chromosomes)[chrom_ind].fill_read(read, 0, read_start, read_chrom_space);

    // Reverse complement if necessary:
    if (reverse) rev_comp(read, read_chrom_space);

    /*
     Adding read with errors:
     */
    uint64 read_pos = 0;
    uint64 current_length = 0;
    uint64 rndi;
    while (current_length < read_length) {
        if (!insertions.empty() && read_pos == insertions.front()) {
            rndi = static_cast<uint64>(runif_01(eng) * 4);
            fastq_pool.push_back(read[read_pos]);
            fastq_pool.push_back(jlp::bases[rndi]);
            insertions.pop_front();
            current_length += 2;
        } else if (!deletions.empty() && read_pos == deletions.front()) {
            deletions.pop_front();
        } else if (!substitutions.empty() && read_pos == substitutions.front()) {
            rndi = static_cast<uint64>(runif_01(eng) * 3);
            fastq_pool.push_back(mm_nucleos[nt_map[read[read_pos]]][rndi]);
            substitutions.pop_front();
            current_length++;
        } else {
            fastq_pool.push_back(read[read_pos]);
            current_length++;
        }
        read_pos++;
    }

    fastq_pool.push_back('\n');
    fastq_pool.push_back('+');
    fastq_pool.push_back('\n');

    // Adding qualities:
    for (uint64 i = 0; i < split_pos; i++) fastq_pool.push_back(qual_left);
    for (uint64 i = split_pos; i < read_length; i++) fastq_pool.push_back(qual_right);
    fastq_pool.push_back('\n');

    return;
}


template <typename T>
template <typename U>
void PacBioOneGenome<T>::append_pool(const std::string& chrom,
                                     U& fastq_pool,
                                     pcg64& eng) {

    // Make sure it has enough memory reserved:
    fastq_pool.reserve(fastq_pool.size() + read_length * 3 + 10);

    // Boolean for whether we take the reverse side:
    bool reverse = runif_01(eng) < 0.5;

    // ID line:
    fastq_pool.push_back('@');
    for (const char& c : name) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    for (const char& c : (*chromosomes)[chrom_ind].name) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    for (const char& c : std::to_string(read_start)) fastq_pool.push_back(c);
    fastq_pool.push_back('-');
    if (reverse) {
        fastq_pool.push_back('R');
    } else fastq_pool.push_back('F');
    fastq_pool.push_back('\n');

    // Fill in read:
    fill_read__(chrom, read, 0, read_start, read_chrom_space);

    // Reverse complement if necessary:
    if (reverse) rev_comp(read, read_chrom_space);

    /*
     Adding read with errors:
     */
    uint64 read_pos = 0;
    uint64 current_length = 0;
    uint64 rndi;
    while (current_length < read_length) {
        if (!insertions.empty() && read_pos == insertions.front()) {
            rndi = static_cast<uint64>(runif_01(eng) * 4);
            fastq_pool.push_back(read[read_pos]);
            fastq_pool.push_back(jlp::bases[rndi]);
            insertions.pop_front();
            current_length += 2;
        } else if (!deletions.empty() && read_pos == deletions.front()) {
            deletions.pop_front();
        } else if (!substitutions.empty() && read_pos == substitutions.front()) {
            rndi = static_cast<uint64>(runif_01(eng) * 3);
            fastq_pool.push_back(mm_nucleos[nt_map[read[read_pos]]][rndi]);
            substitutions.pop_front();
            current_length++;
        } else {
            fastq_pool.push_back(read[read_pos]);
            current_length++;
        }
        read_pos++;
    }

    fastq_pool.push_back('\n');
    fastq_pool.push_back('+');
    fastq_pool.push_back('\n');

    // Adding qualities:
    for (uint64 i = 0; i < split_pos; i++) fastq_pool.push_back(qual_left);
    for (uint64 i = split_pos; i < read_length; i++) fastq_pool.push_back(qual_right);
    fastq_pool.push_back('\n');

    return;
}






// `one_read` method
template <typename U>
void PacBioHaplotypes::one_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng) {


    if (hap == haplotypes->size()) {
        finished = true;
        return;
    }

    if (n_reads_vc[hap][chr] == 0 || hap_chrom_seq.empty()) {

        uint64 new_hap = hap;
        uint64 new_chr = chr;
        for (; new_hap < n_reads_vc.size(); new_hap++) {
            while (n_reads_vc[new_hap][new_chr] == 0) {
                new_chr++;
                if (new_chr == n_reads_vc[new_hap].size()) break;
            }
            if (new_chr < n_reads_vc[new_hap].size()) {
                break;
            } else new_chr = 0;
        }

        hap = new_hap;
        chr = new_chr;

        if (hap == haplotypes->size())  {
            finished = true;
            return;
        }

        hap_chrom_seq = (*haplotypes)[hap][chr].get_chrom_full();
    }

    read_makers[hap].one_read<U>(hap_chrom_seq, chr, fastq_pools, eng);

    n_reads_vc[hap][chr]--;

    return;

}


// `re_read` method (for duplicates)
template <typename U>
void PacBioHaplotypes::re_read(std::vector<U>& fastq_pools, bool& finished, pcg64& eng) {

    if (hap == haplotypes->size()) {
        finished = true;
        return;
    }

    read_makers[hap].re_read<U>(hap_chrom_seq, chr, fastq_pools, eng);

    if (n_reads_vc[hap][chr] > 0) n_reads_vc[hap][chr]--;

    return;

}








/*
 ========================================================================================
 ========================================================================================

 Writing reads

 ========================================================================================
 ========================================================================================
 */




//' PacBio chromosome for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void pacbio_ref_cpp(SEXP ref_genome_ptr,
                      const std::string& out_prefix,
                      const int& compress,
                      const std::string& comp_method,
                      const uint64& n_reads,
                      const uint64& n_threads,
                      const bool& show_progress,
                      const uint64& read_pool_size,
                      const double& prob_dup,
                      const double& scale,
                      const double& sigma,
                      const double& loc,
                      const double& min_read_len,
                      const std::vector<double>& read_probs,
                      const std::vector<uint64>& read_lens,
                      const uint64& max_passes,
                      const std::vector<double>& chi2_params_n,
                      const std::vector<double>& chi2_params_s,
                      const std::vector<double>& sqrt_params,
                      const std::vector<double>& norm_params,
                      const double& prob_thresh,
                      const double& prob_ins,
                      const double& prob_del,
                      const double& prob_subst) {

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

    // For doing multithreaded compression after initial uncompressed run:
    uint64 prog_n = n_reads;
    if (compress > 0 && n_threads > 1) prog_n += (n_reads / 2);
    // Progress bar:
    Progress prog_bar(prog_n, show_progress);

    write_reads_cpp_<PacBioReference>(
        read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size, 1,
        n_threads, compress, comp_method, prog_bar);


    return;
}




//' PacBio chromosome for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void pacbio_hap_cpp(SEXP hap_set_ptr,
                    const std::string& out_prefix,
                    const bool& sep_files,
                    const int& compress,
                    const std::string& comp_method,
                    const uint64& n_reads,
                    const uint64& n_threads,
                    const bool& show_progress,
                    const uint64& read_pool_size,
                    const std::vector<double>& haplotype_probs,
                    const double& prob_dup,
                    const double& scale,
                    const double& sigma,
                    const double& loc,
                    const double& min_read_len,
                    const std::vector<double>& read_probs,
                    const std::vector<uint64>& read_lens,
                    const uint64& max_passes,
                    const std::vector<double>& chi2_params_n,
                    const std::vector<double>& chi2_params_s,
                    const std::vector<double>& sqrt_params,
                    const std::vector<double>& norm_params,
                    const double& prob_thresh,
                    const double& prob_ins,
                    const double& prob_del,
                    const double& prob_subst) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    PacBioHaplotypes read_filler_base;

    if (read_probs.size() == 0) {
        read_filler_base =
            PacBioHaplotypes(*hap_set, haplotype_probs,
                           scale, sigma, loc, min_read_len,
                           max_passes, chi2_params_n, chi2_params_s,
                           sqrt_params, norm_params, prob_thresh,
                           prob_ins, prob_del, prob_subst);
    } else {
        read_filler_base =
            PacBioHaplotypes(*hap_set, haplotype_probs,
                           read_probs, read_lens,
                           max_passes, chi2_params_n, chi2_params_s,
                           sqrt_params, norm_params, prob_thresh,
                           prob_ins, prob_del, prob_subst);
    }

    // For doing multithreaded compression after initial uncompressed run:
    uint64 prog_n = n_reads;
    if (compress > 0 && n_threads > 1) prog_n += (n_reads / 2);
    // Progress bar:
    Progress prog_bar(prog_n, show_progress);

    if (sep_files) {

        write_reads_cpp_sep_files_<PacBioHaplotypes>(
            *hap_set, haplotype_probs,
            read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size, 1,
            n_threads, compress, comp_method, prog_bar);

    } else {

        write_reads_cpp_<PacBioHaplotypes>(
            read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size, 1,
            n_threads, compress, comp_method, prog_bar);

    }


    return;
}



