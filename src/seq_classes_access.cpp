
/*
 ********************************************************

 Retrieving and printing info from Ref* and Var* classes, from R environment.

 ********************************************************
 */


#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque


#include "gemino_types.h"  // integer types
#include "seq_classes_ref.h"  // Ref* classes
#include "seq_classes_var.h"  // Var* classes

using namespace Rcpp;



/*
 ----------------------------------------
 Reference genome methods
 ----------------------------------------
 */



//' Function to print info on a `RefGenome`.
//'
//' Access `RefGenome` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_ref_genome(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    ref_genome->print();
    return;
}


//[[Rcpp::export]]
std::vector<uint32> see_ref_genome_seq_sizes(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<uint32> out = ref_genome->seq_sizes();
    return out;
}

// Check for no duplicates before this fxn
//[[Rcpp::export]]
void remove_ref_genome_seqs(
        SEXP ref_genome_ptr,
        std::vector<uint32> seq_inds) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::deque<RefSequence>& sequences(ref_genome->sequences);

    std::sort(seq_inds.begin(), seq_inds.end());

    for (uint32 i = 0; i < seq_inds.size(); i++) {
        uint32 j = seq_inds[i];
        sequences.erase(sequences.begin() + j);
        for (uint32 k = i+1; k < seq_inds.size(); k++) seq_inds[k]--;
    }
    clear_memory<std::deque<RefSequence>>(sequences);
    return;
}


//[[Rcpp::export]]
std::vector<std::string> see_ref_genome_seq_names(SEXP ref_genome_ptr) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::vector<std::string> out;
    out.reserve(ref_genome->size());
    for (const RefSequence& seq : (*ref_genome).sequences) out.push_back(seq.name);
    return out;
}

// Check for seq_inds and names lengths being equal before this fxn
//[[Rcpp::export]]
void set_ref_genome_seq_names(
        SEXP ref_genome_ptr,
        const std::vector<uint32>& seq_inds,
        const std::vector<std::string>& names) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    for (uint32 i = 0; i < seq_inds.size(); i++) {
        (*ref_genome)[seq_inds[i]].name = names[i];
    }
    return;
}

//[[Rcpp::export]]
std::string see_ref_genome_seq(SEXP ref_genome_ptr, const uint32& seq_ind) {
    XPtr<RefGenome> ref_genome(ref_genome_ptr);
    std::string out = (*ref_genome)[seq_ind].nucleos;
    return out;
}



/*
 ----------------------------------------
 Variants methods
 ----------------------------------------
 */


//' Function to print info on a VarSet.
//'
//' Access `VarSet` class's print method from R.
//'
//' @noRd
//'
//[[Rcpp::export]]
void print_var_set(SEXP var_set_ptr) {
    XPtr<VarSet> var_set(var_set_ptr);
    var_set->print();
    return;
}



