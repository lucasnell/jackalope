
#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng

#include <fstream> // for writing FASTQ files
#include "zlib.h" // for writing to compressed FASTQ
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // integer types
#include "ref_classes.h"  // Ref* classes
#include "hap_classes.h"  // Hap* classes
#include "hts.h"  // generic sequencing classes
#include "hts_illumina.h"  // Illumina-specific classes
#include "str_manip.h"  // rev_comp
#include "io.h"  // File* types
#include "util.h"  // split_int


using namespace Rcpp;





// Sample one set of read strings (each with 4 lines: ID, chromosome, "+", quality)
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::one_read(std::vector<U>& fastq_pools,
                                    bool& finished,
                                    pcg64& eng) {

    /*
     Sample fragment info, and set the chromosome space(s) required for these read(s).
     */
    chrom_indels_frag(eng);

    if (constr_info.chrom_ind == chromosomes->size())  {
        finished = true;
        return;
    }

    // Fill the reads and qualities
    append_pools<U>(fastq_pools, eng);

    return;
}

// Overloaded for when we input a haplotype chromosome stored as string
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::one_read(const std::string& chrom,
                                    const uint64& chrom_i,
                                    std::vector<U>& fastq_pools,
                                    pcg64& eng) {

    constr_info.chrom_ind = chrom_i;

    /*
     Sample indels, and set the chromosome space(s) required for these read(s).
     */
    indels_frag(eng);

    // Fill the reads and qualities
    append_pools<U>(chrom, fastq_pools, eng);

    return;
}

/*
 Same as above, but for a duplicate. It's assumed that `one_read` has been
 run once before.
 */
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::re_read(std::vector<U>& fastq_pools,
                                   bool& finished,
                                   pcg64& eng) {

    // Here I'm just re-doing indels bc it's a duplicate.
    just_indels(eng);

    // Fill the reads and qualities
    append_pools<U>(fastq_pools, eng);

    return;
}
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::re_read(const std::string& chrom,
                                   const uint64& chrom_i,
                                   std::vector<U>& fastq_pools,
                                   pcg64& eng) {

    constr_info.chrom_ind = chrom_i;

    // Here I'm just re-doing indels bc it's a duplicate.
    just_indels(eng);

    // Fill the reads and qualities
    append_pools<U>(chrom, fastq_pools, eng);

    return;
}




// Sample for insertion and deletion positions
template <typename T>
void IlluminaOneGenome<T>::sample_indels(pcg64& eng) {

    const uint64& frag_len(constr_info.frag_len);

    for (uint64 r = 0; r < insertions.size(); r++) {
        uint64 frag_pos = 0;
        uint64 length_now = 0;
        double u;
        std::deque<uint64>& ins(insertions[r]);
        std::deque<uint64>& del(deletions[r]);
        const double& ins_prob(ins_probs[r]);
        const double& del_prob(del_probs[r]);
        ins.clear();
        del.clear();
        while (length_now < read_length && frag_pos < frag_len) {
            u = runif_01(eng);
            if (u > (ins_prob + del_prob)) {
                length_now++;
            } else if (u > ins_prob) {
                del.push_back(frag_pos);
            } else {
                if (length_now == (read_length - 1)) {
                    length_now++;
                } else {
                    ins.push_back(frag_pos);
                    length_now += 2;
                }
            }
            frag_pos++;
        }
    }

    return;
}

// Adjust chromosome spaces
template <typename T>
void IlluminaOneGenome<T>::adjust_chrom_spaces() {

    std::vector<uint64>& read_chrom_spaces(constr_info.read_chrom_spaces);
    std::vector<std::string>& reads(constr_info.reads);
    const uint64& frag_len(constr_info.frag_len);

    for (uint64 r = 0; r < insertions.size(); r++) {
        /*
         I'm adding deletions because more deletions mean that I need
         more chromosome bases to achieve the same read length.
         Insertions means I need fewer.
         */
        sint64 indel_effect = static_cast<sint64>(deletions[r].size()) -
            static_cast<sint64>(insertions[r].size());
        /*
         In addition to indels, below corrects for situation where a small
         fragment size was sampled.
         (Because the indel sampler stops when it reaches the fragment end,
         we don't need to account for that.)
         */
        read_chrom_spaces[r] = std::min(read_length + indel_effect, frag_len);
        // Adjust `reads` so it can hold the necessary chromosome space:
        if (reads[r].size() != read_chrom_spaces[r]) {
            reads[r].resize(read_chrom_spaces[r], 'N');
        }
        // Now including effect of barcode:
        read_chrom_spaces[r] -= constr_info.barcode.size();
    }

    return;
}


/*
 Sample a chromosome, indels, fragment length, and starting position for the fragment.
 Lastly, it sets the chromosome spaces required for these reads.
 */
template <typename T>
void IlluminaOneGenome<T>::chrom_indels_frag(pcg64& eng) {

    uint64& chrom_ind(constr_info.chrom_ind);
    uint64& frag_len(constr_info.frag_len);
    uint64& frag_start(constr_info.frag_start);

    // Get chromosome:
    chrom_ind = 0;
    while (chrom_ind < chrom_reads.size() && chrom_reads[chrom_ind] == 0) chrom_ind++;
    if (chrom_ind == chromosomes->size()) return;

    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample fragment length:
    frag_len = static_cast<uint64>(frag_lengths(eng));
    if (frag_len < frag_len_min) frag_len = frag_len_min;
    if (frag_len > frag_len_max) frag_len = frag_len_max;

    // Sample fragment starting position:
    if (frag_len >= chrom_len) {
        frag_len = chrom_len;
        frag_start = 0;
    } else {
        double u = runif_01(eng);
        frag_start = static_cast<uint64>(u * (chrom_len - frag_len + 1));
    }

    // Sample indels:
    sample_indels(eng);

    // Adjust chromosome spaces:
    adjust_chrom_spaces();

    return;
}



/*
 Sample indels, fragment length, and starting position for the fragment.
 Lastly, it sets the chromosome spaces required for these reads.
 This is for when the chromosome is already set.
 */
template <typename T>
void IlluminaOneGenome<T>::indels_frag(pcg64& eng) {

    uint64& chrom_ind(constr_info.chrom_ind);
    uint64& frag_len(constr_info.frag_len);
    uint64& frag_start(constr_info.frag_start);

    uint64 chrom_len = (*chromosomes)[chrom_ind].size();

    // Sample fragment length:
    frag_len = static_cast<uint64>(frag_lengths(eng));
    if (frag_len < frag_len_min) frag_len = frag_len_min;
    if (frag_len > frag_len_max) frag_len = frag_len_max;

    // Sample fragment starting position:
    if (frag_len >= chrom_len) {
        frag_len = chrom_len;
        frag_start = 0;
    } else {
        double u = runif_01(eng);
        frag_start = static_cast<uint64>(u * (chrom_len - frag_len + 1));
    }

    // Sample indels:
    sample_indels(eng);

    // Adjust chromosome spaces:
    adjust_chrom_spaces();

    return;
}


/*
 Same as `chrom_indels_frag` above, but for duplicates.
 This means skipping the chromosome and fragment info parts.
 */
template <typename T>
void IlluminaOneGenome<T>::just_indels(pcg64& eng) {

    // Sample indels:
    sample_indels(eng);

    // Adjust chromosome spaces:
    adjust_chrom_spaces();

    return;
}


template <typename U>
void fill_fq_lines(U& fq_pool,
                   const std::string& name,
                   const std::string& chrom_name,
                   const std::string& read,
                   const std::string& qual,
                   const uint64& i,
                   const uint64& start,
                   const bool& paired,
                   bool& reverse) {

    // Combine into 4 lines of output per read:
    // ID line:
    fq_pool.push_back('@');
    for (const char& c : name) fq_pool.push_back(c);
    fq_pool.push_back('-');
    for (const char& c : chrom_name) fq_pool.push_back(c);
    fq_pool.push_back('-');
    for (const char& c : std::to_string(start)) fq_pool.push_back(c);
    fq_pool.push_back('-');
    if (reverse) {
        fq_pool.push_back('R');
    } else fq_pool.push_back('F');
    if (paired) {
        fq_pool.push_back('/');
        for (const char& c : std::to_string(i+1)) fq_pool.push_back(c);
    }
    fq_pool.push_back('\n');
    // The rest:
    for (const char& c : read) fq_pool.push_back(c);
    fq_pool.push_back('\n');
    fq_pool.push_back('+');
    fq_pool.push_back('\n');
    for (const char& c : qual) fq_pool.push_back(c);
    fq_pool.push_back('\n');

    // If doing paired reads, the second one should be the reverse of the first
    reverse = !reverse;

    return;

}




/*
 Sample one set of read strings (each with 4 lines: ID, chromosome, "+", quality),
 then append that to the `fastq_pools` vector.
 This function does NOT do anything with fragments.
 That should be done outside this function.
 */
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::append_pools(std::vector<U>& fastq_pools,
                                         pcg64& eng) {

    uint64 n_read_ends = ins_probs.size();
    if (fastq_pools.size() != n_read_ends) fastq_pools.resize(n_read_ends);

    const std::string& barcode(constr_info.barcode);
    const std::vector<uint64>& read_chrom_spaces(constr_info.read_chrom_spaces);

    // Just making this reference to keep lines from getting very long.
    const uint64& chrom_ind(constr_info.chrom_ind);

    // Boolean for whether we take the reverse side first:
    bool reverse = runif_01(eng) < 0.5;
    for (uint64 i = 0; i < n_read_ends; i++) {
        std::string& read(constr_info.reads[i]);
        std::string& qual(constr_info.quals[i]);

        // Read starting location depends on if mate-pair and if reverse strand:
        uint64 start;
        if ((!matepair && !reverse) || (matepair && reverse)) {
            start = constr_info.frag_start;
        } else {
            start = constr_info.frag_start + constr_info.frag_len -
                read_chrom_spaces[i];
        }

        /*
         Filling `read` only depends on whether it's the reverse strand:
        */
        if (!reverse) {
            // Fill in read starting with position after barcode:
            (*(chromosomes))[chrom_ind].fill_read(
                    read, barcode.size(),
                    start, read_chrom_spaces[i]);
        } else {
            /*
             If doing reverse, we can add the actual chromosome to the front of the
             read instead of starting at `barcode.size()`.

             (Even though `rev_comp` requires only T, C, A, G, or N characters,
             `read` should already have filler 'N' chars present at initialization
             of the `constr_info` field, so no need to add any now.)
             */
            (*(chromosomes))[chrom_ind].fill_read(
                    read, 0,
                    start, read_chrom_spaces[i]);
            // Now do reverse complement:
            rev_comp(read);
        }

        // Now fill barcode:
        for (uint64 i = 0; i < barcode.size(); i++) read[i] = barcode[i];

        // Sample mapping quality and add errors to read:
        qual_errors[i].fill_read_qual(read, qual, insertions[i], deletions[i], eng);

        std::string chrom_name = (*chromosomes)[chrom_ind].name;

        // Combine into 4 lines of output per read, and add to `fastq_pools[i]`
        fill_fq_lines<U>(fastq_pools[i], name, chrom_name, read, qual, i, start,
                         paired, reverse);

    }

    if (chrom_reads[constr_info.chrom_ind] < n_read_ends) {
        chrom_reads[constr_info.chrom_ind] = 0;
    } else chrom_reads[constr_info.chrom_ind] -= n_read_ends;

    return;
}

// Overloaded for when a haplotype chromosome sequence is provided
template <typename T>
template <typename U>
void IlluminaOneGenome<T>::append_pools(const std::string& chrom,
                                        std::vector<U>& fastq_pools,
                                        pcg64& eng) {

    uint64 n_read_ends = ins_probs.size();
    if (fastq_pools.size() != n_read_ends) fastq_pools.resize(n_read_ends);

    const std::string& barcode(constr_info.barcode);
    const std::vector<uint64>& read_chrom_spaces(constr_info.read_chrom_spaces);

    // Just making this reference to keep lines from getting very long.
    const uint64& chrom_ind(constr_info.chrom_ind);

    // Boolean for whether we take the reverse side first:
    bool reverse = runif_01(eng) < 0.5;
    for (uint64 i = 0; i < n_read_ends; i++) {
        std::string& read(constr_info.reads[i]);
        std::string& qual(constr_info.quals[i]);

        // Read starting location depends on if mate-pair and if reverse strand:
        uint64 start;
        if ((!matepair && !reverse) || (matepair && reverse)) {
            start = constr_info.frag_start;
        } else {
            start = constr_info.frag_start + constr_info.frag_len -
                read_chrom_spaces[i];
        }

        /*
         Filling `read` only depends on whether it's the reverse strand:
        */
        if (!reverse) {

            // Fill in read starting with position after barcode:
            fill_read__(chrom, read, barcode.size(), start, read_chrom_spaces[i]);
            // `fill_read__` fxn found in hts.h


        } else {
            /*
             If doing reverse, we can add the actual chromosome to the front of the
             read instead of starting at `barcode.size()`.

             (Even though `rev_comp` requires only T, C, A, G, or N characters,
             `read` should already have filler 'N' chars present at initialization
             of the `constr_info` field, so no need to add any now.)
             */
            fill_read__(chrom, read, 0, start, read_chrom_spaces[i]);

            // Now do reverse complement:
            rev_comp(read);
        }

        // Now fill barcode:
        for (uint64 i = 0; i < barcode.size(); i++) read[i] = barcode[i];

        // Sample mapping quality and add errors to read:
        qual_errors[i].fill_read_qual(read, qual, insertions[i], deletions[i], eng);

        std::string chrom_name = (*chromosomes)[chrom_ind].name;

        // Combine into 4 lines of output per read, and add to `fastq_pools[i]`
        fill_fq_lines<U>(fastq_pools[i], name, chrom_name, read, qual, i, start,
                         paired, reverse);

    }

    return;
}











// If only providing rng and id info, sample for a haplotype, then make read(s):
template <typename U>
void IlluminaHaplotypes::one_read(std::vector<U>& fastq_pools,
                                bool& finished,
                                pcg64& eng) {

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
    if (paired && n_reads_vc[hap][chr] > 0) n_reads_vc[hap][chr]--;

    return;
}
/*
 -------------
 `re_read` methods (for duplicates)
 -------------
 */
template <typename U>
void IlluminaHaplotypes::re_read(std::vector<U>& fastq_pools,
                               bool& finished,
                               pcg64& eng) {

    if (hap == haplotypes->size()) {
        finished = true;
        return;
    }

    read_makers[hap].re_read<U>(hap_chrom_seq, chr, fastq_pools, eng);

    if (n_reads_vc[hap][chr] > 0) n_reads_vc[hap][chr]--;
    if (paired && n_reads_vc[hap][chr] > 0) n_reads_vc[hap][chr]--;

    return;
}






/*
 ========================================================================================
 ========================================================================================

 Writing reads

 ========================================================================================
 ========================================================================================
 */









//' Illumina chromosome for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void illumina_ref_cpp(SEXP ref_genome_ptr,
                      const bool& paired,
                      const bool& matepair,
                      const std::string& out_prefix,
                      const int& compress,
                      const std::string& comp_method,
                      const uint64& n_reads,
                      const double& prob_dup,
                      const uint64& n_threads,
                      const bool& show_progress,
                      const uint64& read_pool_size,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint64& frag_len_min,
                      const uint64& frag_len_max,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                      const std::vector<std::vector<std::vector<uint8>>>& quals1,
                      const double& ins_prob1,
                      const double& del_prob1,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                      const std::vector<std::vector<std::vector<uint8>>>& quals2,
                      const double& ins_prob2,
                      const double& del_prob2,
                      const std::vector<std::string>& barcodes) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    IlluminaReference read_filler_base;

    uint64 n_read_ends;
    if (paired) {
        n_read_ends = 2;
        read_filler_base =
            IlluminaReference(*ref_genome, matepair,
                              frag_len_shape, frag_len_scale,
                              frag_len_min, frag_len_max,
                              qual_probs1, quals1, ins_prob1, del_prob1,
                              qual_probs2, quals2, ins_prob2, del_prob2,
                              barcodes[0]);
    } else {
        n_read_ends = 1;
        read_filler_base =
            IlluminaReference(*ref_genome,
                              frag_len_shape, frag_len_scale,
                              frag_len_min, frag_len_max,
                              qual_probs1, quals1, ins_prob1, del_prob1,
                              barcodes[0]);
    }

    // For doing multithreaded compression after initial uncompressed run:
    uint64 prog_n = n_reads;
    if (compress > 0 && n_threads > 1) prog_n += (n_reads / 2);
    // Progress bar:
    Progress prog_bar(prog_n, show_progress);

    write_reads_cpp_<IlluminaReference>(
        read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size,
        n_read_ends, n_threads, compress, comp_method, prog_bar);

    return;
}






//' Illumina chromosome for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void illumina_hap_cpp(SEXP hap_set_ptr,
                      const bool& paired,
                      const bool& matepair,
                      const std::string& out_prefix,
                      const bool& sep_files,
                      const int& compress,
                      const std::string& comp_method,
                      const uint64& n_reads,
                      const double& prob_dup,
                      const uint64& n_threads,
                      const bool& show_progress,
                      const uint64& read_pool_size,
                      const std::vector<double>& haplotype_probs,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint64& frag_len_min,
                      const uint64& frag_len_max,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                      const std::vector<std::vector<std::vector<uint8>>>& quals1,
                      const double& ins_prob1,
                      const double& del_prob1,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                      const std::vector<std::vector<std::vector<uint8>>>& quals2,
                      const double& ins_prob2,
                      const double& del_prob2,
                      const std::vector<std::string>& barcodes) {

    XPtr<HapSet> hap_set(hap_set_ptr);
    IlluminaHaplotypes read_filler_base;

    uint64 n_read_ends;

    if (paired) {

        n_read_ends = 2;
        read_filler_base = IlluminaHaplotypes(*hap_set, haplotype_probs,
                                            matepair,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            qual_probs2, quals2, ins_prob2, del_prob2,
                                            barcodes);

    } else {

        n_read_ends = 1;
        read_filler_base = IlluminaHaplotypes(*hap_set, haplotype_probs,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            barcodes);

    }

    // For doing multithreaded compression after initial uncompressed run:
    uint64 prog_n = n_reads;
    if (compress > 0 && n_threads > 1) prog_n += (n_reads / 2);
    // Progress bar:
    Progress prog_bar(prog_n, show_progress);

    if (sep_files) {

        write_reads_cpp_sep_files_<IlluminaHaplotypes>(
            *hap_set, haplotype_probs,
            read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size,
            n_read_ends, n_threads, compress, comp_method, prog_bar);


    } else {

        write_reads_cpp_<IlluminaHaplotypes>(
            read_filler_base, out_prefix, n_reads, prob_dup, read_pool_size,
            n_read_ends, n_threads, compress, comp_method, prog_bar);

    }

    return;
}






