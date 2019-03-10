
#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng

#include <fstream> // for writing FASTQ files
#include <zlib.h>  // for writing to compressed FASTQ

#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "sequencer.h"  // SequenceIdentifierInfo class
#include "illumina.h"  // Illumina-specific classes

using namespace Rcpp;




//' Illumina sequence for reference object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
void illumina_ref_cpp(SEXP ref_genome_ptr,
                      const std::string& out_prefix,
                      const bool& compress,
                      const uint32& n_reads,
                      const double& prob_pcr_dup,
                      const uint32& n_cores,
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
                      const std::string& barcode,
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

    if (qual_probs2.size() == 0) {
        read_filler_base =
            IlluminaReference(*ref_genome, frag_len_shape, frag_len_scale,
                              frag_len_min, frag_len_max,
                              qual_probs1, quals1, ins_prob1, del_prob1,
                              barcode);

    } else {
        read_filler_base =
            IlluminaReference(*ref_genome, frag_len_shape, frag_len_scale,
                              frag_len_min, frag_len_max,
                              qual_probs1, quals1, ins_prob1, del_prob1,
                              qual_probs2, quals2, ins_prob2, del_prob2,
                              barcode);
    }

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);

    if (compress) {
        illumina_cpp_<IlluminaReference, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_pcr_dup,
                read_chunk_size, n_cores);
    } else {
        illumina_cpp_<IlluminaReference, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_pcr_dup,
                read_chunk_size, n_cores);
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
                      const std::string& out_prefix,
                      const bool& compress,
                      const uint32& n_reads,
                      const double& prob_pcr_dup,
                      const uint32& n_cores,
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

    if (qual_probs2.size() == 0) {
        read_filler_base = IlluminaVariants(*var_set, variant_probs,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            barcodes);
    } else {
        read_filler_base = IlluminaVariants(*var_set, variant_probs,
                                            frag_len_shape, frag_len_scale,
                                            frag_len_min, frag_len_max,
                                            qual_probs1, quals1, ins_prob1, del_prob1,
                                            qual_probs2, quals2, ins_prob2, del_prob2,
                                            barcodes);
    }

    SequenceIdentifierInfo ID_info_base(instrument, run_number, flowcell_ID, lane,
                                        tile, x_pos, y_pos, read, is_filtered,
                                        control_number, sample_number);


    if (compress) {
        write_illumina_<IlluminaVariants, gzFile>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_pcr_dup,
                read_chunk_size, n_cores);
    }else {
        write_illumina_<IlluminaVariants, std::ofstream>(
                read_filler_base, ID_info_base, out_prefix, n_reads, prob_pcr_dup,
                read_chunk_size, n_cores);
    }

    return;
}
