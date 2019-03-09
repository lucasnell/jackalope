
#include <RcppArmadillo.h>
#include <algorithm> // lower_bound
#include <vector>  // vector class
#include <string>  // string class
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // distributions

#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes
#include "sequencer.h"  // SequenceIdentifierInfo class
#include "illumina.h"  // Illumina-specific classes

using namespace Rcpp;

//' Create an XPtr to object that creates paired-end Illumina reads of a reference object.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP create_ref_ill_pe(
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
        const double& del_prob2) {
    XPtr<ReferenceIllumina> xptr (new ReferenceIllumina(frag_len_shape, frag_len_scale,
                                                        frag_len_min, frag_len_max,
                                                        qual_probs1, quals1,
                                                        ins_prob1, del_prob1,
                                                        qual_probs2, quals2,
                                                        ins_prob2, del_prob2));
    return xptr;
}


//' Create an XPtr to object that creates single-end Illumina reads of a reference object.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP create_ref_ill_se(
        const double& frag_len_shape,
        const double& frag_len_scale,
        const uint32& frag_len_min,
        const uint32& frag_len_max,
        const std::vector<std::vector<std::vector<double>>>& qual_probs,
        const std::vector<std::vector<std::vector<uint8>>>& quals,
        const double& ins_prob,
        const double& del_prob) {
    XPtr<ReferenceIllumina> xptr(new ReferenceIllumina(frag_len_shape, frag_len_scale,
                                                       frag_len_min, frag_len_max,
                                                       qual_probs, quals,
                                                       ins_prob, del_prob));
    return xptr;
}


//' Create an XPtr to object that creates paired-end Illumina reads of a variants object.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP create_var_ill_pe(
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
        const double& del_prob2) {
    XPtr<VariantIllumina> xptr(new VariantIllumina(variant_probs,
                                                   frag_len_shape, frag_len_scale,
                                                   frag_len_min, frag_len_max,
                                                   qual_probs1, quals1,
                                                   ins_prob1, del_prob1,
                                                   qual_probs2, quals2,
                                                   ins_prob2, del_prob2));
    return xptr;
}




//' Create an XPtr to object that creates single-end Illumina reads of a variants object.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP create_var_ill_se(
        const std::vector<double>& variant_probs,
        const double& frag_len_shape,
        const double& frag_len_scale,
        const uint32& frag_len_min,
        const uint32& frag_len_max,
        const std::vector<std::vector<std::vector<double>>>& qual_probs,
        const std::vector<std::vector<std::vector<uint8>>>& quals,
        const double& ins_prob,
        const double& del_prob) {
    XPtr<VariantIllumina> xptr(new VariantIllumina(variant_probs,
                                                   frag_len_shape, frag_len_scale,
                                                   frag_len_min, frag_len_max,
                                                   qual_probs, quals,
                                                   ins_prob, del_prob));
    return xptr;
}



