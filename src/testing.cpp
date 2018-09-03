
/*
 ********************************************************

 Methods that are only used for testing.

 ********************************************************
 */


#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <deque>  // deque


#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes

using namespace Rcpp;



//' Make a RefGenome object from a set of sequences
//'
//' @noRd
//[[Rcpp::export]]
SEXP make_ref_genome(const std::deque<std::string>& seqs) {
    XPtr<RefGenome> ref_genome(new RefGenome(seqs), true);
    return ref_genome;
}


//' Make a VarSet object from a RefGenome pointer and # variants
//'
//' @noRd
//[[Rcpp::export]]
SEXP make_var_set(SEXP ref_genome_ptr, const uint32& n_vars) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    XPtr<VarSet> var_set(new VarSet(*ref_genome, n_vars), true);
    return var_set;
}

//' Function to piece together the strings for all sequences in a VarGenome.
//'
//' @noRd
//[[Rcpp::export]]
std::vector<std::string> see_var_genome(SEXP var_set_ptr, const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
    const VarGenome& var_genome((*var_set)[var_ind]);

    std::vector<std::string> out(var_genome.size(), "");
    for (uint32 i = 0; i < var_genome.size(); i++) {
        const VarSequence& var_seq(var_genome[i]);
        out[i] = var_seq.get_seq_full();
    }
    return out;
}

//' See all sequence sizes in a VarSet object.
//'
//' @noRd
//[[Rcpp::export]]
std::vector<uint32> see_sizes(SEXP var_set_ptr, const uint32& var_ind) {

    XPtr<VarSet> var_set(var_set_ptr);
    const VarGenome& var_genome((*var_set)[var_ind]);

    std::vector<uint32> out(var_genome.size());
    for (uint32 i = 0; i < var_genome.size(); i++) {
        const VarSequence& var_seq(var_genome.var_genome[i]);
        out[i] = var_seq.seq_size;
    }
    return out;
}









