//
// Alter reference genome sequences
//



#include <RcppArmadillo.h>

#include <algorithm>  // random_shuffle
#include <deque>
#include <string>
#include <vector>

#include "sequence_classes.h"   // new classes
#include "gemino_types.h"  // integer types


using namespace Rcpp;


namespace alter_scaffs {
    // wrapper around R's RNG such that we get a uniform distribution over
    // [0,n) as required by the STL algorithm
    // (see http://gallery.rcpp.org/articles/stl-random-shuffle/)
    inline int rand_wrapper(const int n) { return std::floor(unif_rand()*n); }
}



// ======================================================================================
// ======================================================================================

//  Merge sequences

// ======================================================================================
// ======================================================================================



//' Merge a reference genome into a single sequence.
//'
//'
//' @param ref_ An external pointer (R class \code{externalptr}) to a
//'     \code{RefGenome} class in C++ (the full class in C++ is
//'     \code{Rcpp::XPtr<RefGenome>}).
//'
//' @return Nothing. Changes are made in place.
//'
//' @name merge_sequences
//'
//' @noRd
//'
//[[Rcpp::export]]
void merge_sequences(SEXP ref_) {

    XPtr<RefGenome> reference(ref_);
    std::deque<RefSequence>& seqs(reference->sequences);

    // Shuffling reference info.
    std::random_shuffle(seqs.begin(), seqs.end(), alter_scaffs::rand_wrapper);

    // Merging the back sequences to the first one:
    std::string& nts(seqs.front().nucleos);
    reference->old_names.push_back(seqs.front().name);
    seqs.front().name = "MERGE";
    uint32 i = seqs.size() - 1;
    while (seqs.size() > 1) {
        nts += seqs[i].nucleos;
        reference->old_names.push_back(seqs[i].name);
        --i;
        seqs.pop_back();
    }
    // clear memory in string
    std::string(nts.begin(), nts.end()).swap(nts);
    // clear memory in deque
    std::deque<RefSequence>(seqs.begin(), seqs.end()).swap(seqs);

    reference->merged = true;

    return;
}









// ======================================================================================
// ======================================================================================

//  Filter sequences

// ======================================================================================
// ======================================================================================

//' Filter reference genome sequences by size or for a proportion of total nucleotides.
//'
//'
//' @param ref_ An external pointer (R class \code{externalptr}) to a
//'     \code{RefGenome} class in C++ (the full class in C++ is
//'     \code{Rcpp::XPtr<RefGenome>}).
//' @param min_seq_size Integer minimum sequence size to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//' @param out_seq_prop Numeric proportion of total sequence to keep.
//'     Defaults to \code{0}, which results in this argument being ignored.
//'
//' @return Nothing. Changes are made in place.
//'
//' @name filter_sequences
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
void filter_sequences(SEXP ref_,
                      const uint32& min_seq_size = 0,
                      const double& out_seq_prop = 0) {

    XPtr<RefGenome> reference(ref_);
    std::deque<RefSequence>& seqs(reference->sequences);

    // Checking for sensible inputs
    if (out_seq_prop <= 0 && min_seq_size == 0) {
        stop("Specify > 0 for min_seq_size or out_seq_prop");
    }
    if (out_seq_prop > 0 && min_seq_size > 0) {
        stop("Specify > 0 for min_seq_size OR out_seq_prop");
    }
    if (out_seq_prop > 1) stop("out_seq_prop must be between 0 and 1");

    // Sorting sequence set by size (largest first)
    std::sort(seqs.begin(), seqs.end(), std::greater<RefSequence>());

    // Index that will point to the first sequence to be deleted
    uint32 i = 0;
    // Keeping track of total genome size after filtering
    double out_seq = 0;

    if (min_seq_size > 0) {
        if (seqs.back().size() >= min_seq_size) return;
        if (seqs[i].size() < min_seq_size) {
            stop("Desired minimum scaffold size is too large. None found. "
                     "The minimum size is " + std::to_string(seqs[i].size())
            );
        }
        // after below, `iter` points to the first sequence smaller than the minimum
        while (seqs[i].size() >= min_seq_size) {
            out_seq += static_cast<double>(seqs[i].size());
            ++i;
        }
    } else {
        // Changing total_size to double so I don't have to worry about integer division
        // being a problem
        double total_seq = static_cast<double>(reference->total_size);
        out_seq = static_cast<double>(seqs[i].size());
        while (out_seq / total_seq < out_seq_prop) {
            ++i;
            out_seq += static_cast<double>(seqs[i].size());
        }
        // Getting `i` to point to the first item to be deleted:
        ++i;
    }

    // Erasing using `iter`
    if (i < seqs.size()) {
        seqs.erase(seqs.begin() + i, seqs.end());
        // clear memory:
        std::deque<RefSequence>(seqs.begin(), seqs.end()).swap(seqs);
    }

    reference->total_size = static_cast<uint64>(out_seq);

    return;
}
