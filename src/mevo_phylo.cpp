
/*
 ********************************************************

 Methods for molecular evolution using phylogenies

 ********************************************************
 */



#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <progress.hpp>  // for the progress bar


#include "gemino_types.h"  // integer types
#include "sequence_classes.h"  // Var* and Ref* classes
#include "mevo.h"  // samplers
#include "table_sampler.h" // table sampling
#include "pcg.h" // pcg sampler types
#include "mevo_rate_matrices.h"
#include "mevo_phylo.h"


using namespace Rcpp;



// T must be ChunkPhyloSeq or PhyloSeq

template <typename T>
void fill_one_seq(const List& genome_phylo_info, const uint32& i,
                  std::vector<T>& all_seqs) {

    std::string err_msg;

    const List& seq_phylo_info(genome_phylo_info[i]);

    uint32 n_trees = seq_phylo_info.size();
    if (n_trees == 0) {
        err_msg = "\nNo trees supplied on sequence " + std::to_string(i+1);
        throw(Rcpp::exception(err_msg.c_str(), false));
    }

    std::vector<uint32> n_bases_(n_trees);
    std::vector<std::vector<double>> branch_lens_(n_trees);
    std::vector<arma::Mat<uint32>> edges_(n_trees);
    std::vector<std::vector<std::string>> tip_labels_(n_trees);

    for (uint32 j = 0; j < n_trees; j++) {

        const List& phylo_info(seq_phylo_info[j]);
        std::vector<double> branch_lens = as<std::vector<double>>(
            phylo_info["branch_lens"]);
        arma::Mat<uint32> edges = as<arma::Mat<uint32>>(phylo_info["edges"]);
        std::vector<std::string> tip_labels = as<std::vector<std::string>>(
            phylo_info["labels"]);
        if (branch_lens.size() != edges.n_rows) {
            err_msg = "\nBranch lengths and edges don't have the same ";
            err_msg += "size on sequence " + std::to_string(i+1) + " and tree ";
            err_msg += std::to_string(j+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }
        if (branch_lens.size() == 0) {
            err_msg = "\nEmpty tree on sequence " + std::to_string(i+1);
            err_msg += " and tree " + std::to_string(j+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }
        uint32 start = as<uint32>(phylo_info["start"]);
        sint64 end = as<sint64>(phylo_info["end"]);
        if (end < start) {
            err_msg = "\nEnd position < start position on sequence ";
            err_msg += std::to_string(i+1) + " and tree " + std::to_string(j+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        n_bases_[j] = end - start + 1;
        branch_lens_[j] = branch_lens;
        edges_[j] = edges;
        tip_labels_[j] = tip_labels;
    }


    all_seqs[i] = T(n_bases_, branch_lens_, edges_, tip_labels_);

    return;

}




//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees(List genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<std::vector<PhyloSeq>> all_seqs_xptr(
            new std::vector<PhyloSeq>(n_seqs)
    );
    std::vector<PhyloSeq>& all_seqs(*all_seqs_xptr);

    for (uint32 i = 0; i < n_seqs; i++) {
        fill_one_seq<PhyloSeq>(genome_phylo_info, i, all_seqs);
    }

    return all_seqs_xptr;
}


//' Create XPtr to nested vector of PhyloTree objects from phylogeny information.
//'
//' Same as above, but chunked.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP phylo_info_to_trees_chunk(List genome_phylo_info) {

    uint32 n_seqs = genome_phylo_info.size();

    if (n_seqs == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    XPtr<std::vector<ChunkPhyloSeq>> all_seqs_xptr(
            new std::vector<ChunkPhyloSeq>(n_seqs)
    );
    std::vector<ChunkPhyloSeq>& all_seqs(*all_seqs_xptr);

    for (uint32 i = 0; i < n_seqs; i++) {
        fill_one_seq<ChunkPhyloSeq>(genome_phylo_info, i, all_seqs);
    }

    return all_seqs_xptr;
}








