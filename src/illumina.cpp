
#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng

#include <fstream> // for writing FASTQ files
#include <zlib.h>  // for writing to compressed FASTQ
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "gemino_types.hpp"  // integer types
#include "seq_classes_ref.hpp"  // Ref* classes
#include "seq_classes_var.hpp"  // Var* classes
#include "sequencer.hpp"  // SequenceIdentifierInfo class
#include "illumina.hpp"  // Illumina-specific classes

using namespace Rcpp;





// template <typename T>
// IlluminaOneGenome<T>::

// Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality)
template <typename T>
void IlluminaOneGenome<T>::one_read(std::vector<std::string>& fastq_chunks,
                                    pcg64& eng,
                                    SequenceIdentifierInfo& ID_info) {

    /*
     Sample fragment info, and set the sequence space(s) required for these read(s).
     */
    seq_indels_frag(eng);

    // Fill the reads and qualities
    append_chunks(fastq_chunks, eng, ID_info);

    return;
}

/*
 Same as above, but for a duplicate. It's assumed that `one_read` has been
 run once before.
 */
template <typename T>
void IlluminaOneGenome<T>::re_read(std::vector<std::string>& fastq_chunks,
                                   pcg64& eng,
                                   SequenceIdentifierInfo& ID_info) {

    // Here I'm just re-doing indels bc it's a duplicate.
    just_indels(eng);

    // Fill the reads and qualities
    append_chunks(fastq_chunks, eng, ID_info);

    return;
}

/*
 Add information about a RefGenome or VarGenome object
 This is used when making multiple samplers that share most info except for
 that related to the sequence object.
 */
template <typename T>
void IlluminaOneGenome<T>::add_seq_info(const T& seq_object, const std::string& barcode) {

    seq_lengths = seq_object.seq_sizes();
    sequences = &seq_object;

    std::vector<double> probs_;
    probs_.reserve(seq_lengths.size());
    for (uint i = 0; i < seq_lengths.size(); i++) {
        probs_.push_back(static_cast<double>(seq_lengths[i]));
    }
    seq_sampler = AliasSampler(probs_);

    constr_info.barcode = barcode;

    return;
}


// Construct sequence-sampling probabilities:
template <typename T>
void IlluminaOneGenome<T>::construct_seqs() {
    std::vector<double> probs_;
    probs_.reserve(seq_lengths.size());
    for (uint i = 0; i < seq_lengths.size(); i++) {
        probs_.push_back(static_cast<double>(seq_lengths[i]));
    }
    seq_sampler = AliasSampler(probs_);
    return;
}


// Sample for insertion and deletion positions
template <typename T>
void IlluminaOneGenome<T>::sample_indels(pcg64& eng) {

    const uint32& frag_len(constr_info.frag_len);

    for (uint32 r = 0; r < insertions.size(); r++) {
        uint32 frag_pos = 0;
        uint32 length_now = 0;
        double u;
        std::deque<uint32>& ins(insertions[r]);
        std::deque<uint32>& del(deletions[r]);
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

// Adjust sequence spaces
template <typename T>
void IlluminaOneGenome<T>::adjust_seq_spaces() {

    std::vector<uint32>& read_seq_spaces(constr_info.read_seq_spaces);
    std::vector<std::string>& reads(constr_info.reads);
    const uint32& frag_len(constr_info.frag_len);

    for (uint32 r = 0; r < insertions.size(); r++) {
        /*
         I'm adding deletions because more deletions mean that I need
         more sequence bases to achieve the same read length.
         Insertions means I need fewer.
         */
        sint32 indel_effect = deletions[r].size() - insertions[r].size();
        /*
         In addition to indels, below corrects for situation where a small
         fragment size was sampled.
         (Because the indel sampler stops when it reaches the fragment end,
         we don't need to account for that.)
         */
        read_seq_spaces[r] = std::min(read_length + indel_effect, frag_len);
        // Adjust `reads` so it can hold the necessary sequence space:
        if (reads[r].size() != read_seq_spaces[r]) {
            reads[r].resize(read_seq_spaces[r], 'N');
        }
        // Now including effect of barcode:
        read_seq_spaces[r] -= constr_info.barcode.size();
    }

    return;
}


/*
 Sample a sequence, indels, fragment length, and starting position for the fragment.
 Lastly, it sets the sequence spaces required for these reads.
 */
template <typename T>
void IlluminaOneGenome<T>::seq_indels_frag(pcg64& eng) {

    uint32& seq_ind(constr_info.seq_ind);
    uint32& frag_len(constr_info.frag_len);
    uint32& frag_start(constr_info.frag_start);

    // Sample sequence:
    seq_ind = seq_sampler.sample(eng);
    uint32 seq_len = (*sequences)[seq_ind].size();

    // Sample fragment length:
    frag_len = static_cast<uint32>(frag_lengths(eng));
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
    sample_indels(eng);

    // Adjust sequence spaces:
    adjust_seq_spaces();

    return;
}


/*
 Same as above, but for duplicates.
 This means skipping the sequence and fragment info parts.
 */
template <typename T>
void IlluminaOneGenome<T>::just_indels(pcg64& eng) {

    // Sample indels:
    sample_indels(eng);

    // Adjust sequence spaces:
    adjust_seq_spaces();

    return;
}



/*
 Sample one set of read strings (each with 4 lines: ID, sequence, "+", quality),
 then append that to the `fastq_chunks` vector.
 This function does NOT do anything with fragments.
 That should be done outside this function.
 */
template <typename T>
void IlluminaOneGenome<T>::append_chunks(std::vector<std::string>& fastq_chunks,
                                         pcg64& eng,
                                         SequenceIdentifierInfo& ID_info) {

    uint32 n_read_ends = ins_probs.size();
    if (fastq_chunks.size() != n_read_ends) fastq_chunks.resize(n_read_ends);

    const std::string& barcode(constr_info.barcode);
    const std::vector<uint32>& read_seq_spaces(constr_info.read_seq_spaces);

    // Boolean for whether we take the reverse side first:
    bool reverse = runif_01(eng) < 0.5;
    for (uint32 i = 0; i < n_read_ends; i++) {
        ID_info.read = i + 1;
        std::string& read(constr_info.reads[i]);
        std::string& qual(constr_info.quals[i]);

        // Read starting location depends on if mate-pair and if reverse strand:
        uint32 start;
        if ((!matepair && !reverse) || (matepair && reverse)) {
            start = constr_info.frag_start;
        } else {
            start = constr_info.frag_start + constr_info.frag_len -
                read_seq_spaces[i];
        }

        /*
         Filling `read` only depends on whether it's the reverse strand:
        */
        if (!reverse) {
            // Fill in read starting with position after barcode:
            (*(sequences))[constr_info.seq_ind].fill_read(
                    read, barcode.size(),
                    start, read_seq_spaces[i]);
        } else {
            /*
             If doing reverse, we can add the actual sequence to the front of the
             read instead of starting at `barcode.size()`.

             (Even though `rev_comp` requires only T, C, A, G, or N characters,
             `read` should already have filler 'N' chars present at initialization
             of the `constr_info` field, so no need to add any now.)
             */
            (*(sequences))[constr_info.seq_ind].fill_read(
                    read, 0,
                    start, read_seq_spaces[i]);
            // Now do reverse complement:
            rev_comp(read);
        }

        // Now fill barcode:
        for (uint i = 0; i < barcode.size(); i++) read[i] = barcode[i];

        // Sample mapping quality and add errors to read:
        qual_errors[i].fill_read_qual(read, qual, insertions[i], deletions[i], eng);

        // If doing paired reads, the second one should be the reverse of the first
        reverse = !reverse;

        // Combine into 4 lines of output per read:
        fastq_chunks[i] += ID_info.get_line() + '\n' + read + "\n+\n" + qual + '\n';
    }

    return;
}







/*
 ========================================================================================
 ========================================================================================

 Writing reads

 ========================================================================================
 ========================================================================================
 */









//' Illumina sequence for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void illumina_ref_cpp(SEXP ref_genome_ptr,
                      const bool& paired,
                      const bool& matepair,
                      const std::string& out_prefix,
                      const bool& compress,
                      const uint32& n_reads,
                      const double& prob_dup,
                      const uint32& n_cores,
                      const bool& show_progress,
                      const uint32& read_chunk_size,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint32& frag_len_min,
                      const uint32& frag_len_max,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                      const std::vector<std::vector<std::vector<uint8>>>& quals1,
                      const double& ins_prob1,
                      const double& del_prob1,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                      const std::vector<std::vector<std::vector<uint8>>>& quals2,
                      const double& ins_prob2,
                      const double& del_prob2,
                      const std::vector<std::string>& barcodes,
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
    IlluminaReference read_filler_base;

    uint32 n_read_ends;
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

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);

    if (compress) {
        write_reads_cpp_<IlluminaReference, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_dup,
                read_chunk_size, n_read_ends, n_cores, show_progress);
    } else {
        write_reads_cpp_<IlluminaReference, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_dup,
                read_chunk_size, n_read_ends, n_cores, show_progress);
    }


    return;
}






//' Illumina sequence for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void illumina_var_cpp(SEXP var_set_ptr,
                      const bool& paired,
                      const bool& matepair,
                      const std::string& out_prefix,
                      const bool& compress,
                      const uint32& n_reads,
                      const double& prob_dup,
                      const uint32& n_cores,
                      const bool& show_progress,
                      const uint32& read_chunk_size,
                      const std::vector<double>& variant_probs,
                      const double& frag_len_shape,
                      const double& frag_len_scale,
                      const uint32& frag_len_min,
                      const uint32& frag_len_max,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs1,
                      const std::vector<std::vector<std::vector<uint8>>>& quals1,
                      const double& ins_prob1,
                      const double& del_prob1,
                      const std::vector<std::vector<std::vector<double>>>& qual_probs2,
                      const std::vector<std::vector<std::vector<uint8>>>& quals2,
                      const double& ins_prob2,
                      const double& del_prob2,
                      const std::vector<std::string>& barcodes,
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
    IlluminaVariants read_filler_base;

    uint32 n_read_ends;
    if (paired) {
        n_read_ends = 2;
        read_filler_base = IlluminaVariants(*var_set, variant_probs,
                                            matepair,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            qual_probs2, quals2, ins_prob2, del_prob2,
                                            barcodes);
    } else {
        n_read_ends = 1;
        read_filler_base = IlluminaVariants(*var_set, variant_probs,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            barcodes);
    }

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);


    if (compress) {
        write_reads_cpp_<IlluminaVariants, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_dup,
                read_chunk_size, n_read_ends, n_cores, show_progress);
    }else {
        write_reads_cpp_<IlluminaVariants, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_dup,
                read_chunk_size, n_read_ends, n_cores, show_progress);
    }

    return;
}
