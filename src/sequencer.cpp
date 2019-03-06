

#include "sequencer.h"  // sequencer classes
#include "gemino_types.h"  // uint32
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes




//' Make pointer to object that simulates Illumina sequencing of a reference genome.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_illumina_ref(SEXP ref_genome_ptr,
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
                       const uint32& read_length,
                       const bool& paired) {

    XPtr<RefGenome> ref_xptr(ref_genome_ptr);
    RefGenome& seq_object(*ref_xptr);

    XPtr<ReferenceIlluminaWGS> illumina_ref(
            new ReferenceIlluminaWGS(seq_object, frag_len_probs, frag_len_region_len,
                                     mis_probs, ins_probs, del_probs, error_region_len,
                                     ins_length_probs, del_length_probs,
                                     qual_means, qual_sds, qual_region_len,
                                     mis_qual_means, mis_qual_sds, mis_qual_region_len,
                                     read_length, paired));

    return illumina_ref;
}


//' Make pointer to object that simulates Illumina sequencing of a variant set object.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_illumina_var(SEXP var_set_ptr,
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
                       const uint32& read_length,
                       const bool& paired) {

    XPtr<VarSet> var_set_xptr(var_set_ptr);
    VarSet& var_set(*var_set_xptr);

    if (variant_probs.size() != var_set.size()) {
        stop("variant_probs.size() != var_set.size()");
    }

    XPtr<VariantIlluminaWGS> illumina_var(
            new VariantIlluminaWGS(var_set, variant_probs,
                                   frag_len_probs, frag_len_region_len,
                                   mis_probs, ins_probs, del_probs, error_region_len,
                                   ins_length_probs, del_length_probs,
                                   qual_means, qual_sds, qual_region_len,
                                   mis_qual_means, mis_qual_sds, mis_qual_region_len,
                                   read_length, paired));

    return illumina_var;
}





//' Make pointer to object that simulates long-read sequencing of a reference genome.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_longread_ref(SEXP ref_genome_ptr,
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
                       const uint32& mis_qual_region_len) {

    XPtr<RefGenome> ref_xptr(ref_genome_ptr);
    RefGenome& seq_object(*ref_xptr);

    XPtr<ReferenceLongReadWGS> longread_ref(
            new ReferenceLongReadWGS(seq_object, frag_len_probs, frag_len_region_len,
                                     mis_probs, ins_probs, del_probs, error_region_len,
                                     ins_length_probs, del_length_probs,
                                     qual_means, qual_sds, qual_region_len,
                                     mis_qual_means, mis_qual_sds, mis_qual_region_len));

    return longread_ref;
}



//' Make pointer to object that simulates long-read sequencing of a reference genome.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP make_longread_var(SEXP var_set_ptr,
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
                       const uint32& mis_qual_region_len) {

    XPtr<VarSet> var_set_xptr(var_set_ptr);
    VarSet& var_set(*var_set_xptr);

    if (variant_probs.size() != var_set.size()) {
        stop("variant_probs.size() != var_set.size()");
    }

    XPtr<VariantLongReadWGS> longread_var(
            new VariantLongReadWGS(var_set, variant_probs,
                                   frag_len_probs, frag_len_region_len,
                                   mis_probs, ins_probs, del_probs, error_region_len,
                                   ins_length_probs, del_length_probs,
                                   qual_means, qual_sds, qual_region_len,
                                   mis_qual_means, mis_qual_sds, mis_qual_region_len));

    return longread_var;
}
