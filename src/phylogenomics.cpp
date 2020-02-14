
/*
 ********************************************************

 Methods for evolving chromosomes along phylogenies / gene trees

 ********************************************************
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <vector>  // vector class
#include <string>  // string class
#include <algorithm>  // lower_bound, sort
#include <deque>  // deque
#include <progress.hpp>  // for the progress bar
#ifdef _OPENMP
#include <omp.h>  // omp
#endif

#include "jackalope_types.h"  // integer types
#include "hap_classes.h"  // Hap* classes
#include "mutator.h"  // TreeMutator
#include "pcg.h" // pcg sampler types
#include "phylogenomics.h"
#include "util.h"  // thread_check


using namespace Rcpp;







/*
 Process one phylogenetic tree for a single chromosome with no recombination.

 Note that this function should be changed if any of these HapChroms differ from
 each other (within the range specified if recombination = true).
 They can already have mutations, but to start out, they must all be the same.
 */
int PhyloOneChrom::one_tree(const uint64& idx,
                            pcg64& eng,
                            Progress& prog_bar) {


    PhyloTree& tree(trees[idx]);

    // For when/if user interrupts function:
    int status;

    // Reset rates for tips:
    status = reset(tree, eng, prog_bar);
    if (status < 0) return status;


    uint64 b1, b2;
    uint64 mut_i = 0;
    sint64 size_mod;
    double b_len;

    /*
     Now iterate through the phylogeny:
     */
    for (uint64 i = 0; i < tree.n_edges; i++) {

        // Checking for abort every edge:
        if (prog_bar.is_aborted() || prog_bar.check_abort()) return -1;

        // Indices for nodes/tips that the branch length in `branch_lens` refers to
        b1 = tree.edges(i,0);
        b2 = tree.edges(i,1);

        // HapChrom we're changing:
        HapChrom& chrom2(*(tip_chroms[b2]));

        if (b1 != b2) {
            // HapChrom object that parent node refers to:
            HapChrom& chrom1(*(tip_chroms[b1]));

            // Update rate indices:
            rates[b2] = rates[b1];

            /*
             Update HapChrom objects for this branch.
             */
            if (idx > 0) mut_i = trees[idx - 1].mut_ends[b1];
            size_mod = chrom2.add_to_back(chrom1, mut_i);

            // Update end point for these new mutations:
            tree.ends[b2] += size_mod;

        }


        // This happens if it's been totally deleted:
        if (tree.starts[b2] == tree.ends[b2]) continue;

        /*
         Now mutate along branch length:
         */
        b_len = tree.branch_lens[i];
#ifdef __JACKALOPE_DIAGNOSTICS
        Rcout << "** b_len " << b_len << std::endl;
#endif
        status = mutator.mutate(b_len, chrom2, eng, prog_bar,
                                tree.starts[b2], tree.ends[b2], rates[b2]);
        if (status < 0) return status;

    }

    // Update indices (non-inclusive) for end of this tree's mutations in
    // `HapChrom::mutations`:
    for (uint64 i = 0; i < tree.n_tips; i++) {
        tree.mut_ends[i] = tip_chroms[i]->mutations.size();
    }

    // Update progress bar:
    prog_bar.increment(tree.end - tree.start + 1);

    return 0;

}







/*
 Reset for a new tree:
 */
int PhyloOneChrom::reset(const PhyloTree& tree,
                         pcg64& eng,
                         Progress& prog_bar) {

    const uint64& start(tree.start);
    const uint64& end(tree.end);

    if (tree.n_tips == 0) {
        throw(Rcpp::exception("\n# tips == zero is non-sensical.", false));
    }

    // Create rates:
    if (rates.size() != n_tips) rates.resize(n_tips);
    uint64 root = tree.edges(0,0); // <-- should be index to root of tree
    // Generate rates for root of tree (`status` is -1 if user interrupts process):
    int status = mutator.new_rates(start, end, rates[root], eng, prog_bar);
    // The rest of the nodes/tips will have rates based on parent nodes
    // as we progress through the tree.

    return status;
}








void PhyloOneChrom::fill_tree_mutator(const List& genome_phylo_info,
                                      const uint64& i,
                                      const TreeMutator& mutator_base) {

    std::string err_msg;

    const List& chrom_phylo_info(genome_phylo_info[i]);

    uint64 n_trees = chrom_phylo_info.size();
    if (n_trees == 0) {
        err_msg = "\nNo trees supplied on chromosome " + std::to_string(i+1);
        throw(Rcpp::exception(err_msg.c_str(), false));
    }

    std::vector<uint64> n_bases_(n_trees);
    std::vector<std::vector<double>> branch_lens_(n_trees);
    std::vector<arma::Mat<uint64>> edges_(n_trees);
    std::vector<std::vector<std::string>> tip_labels_(n_trees);

    for (uint64 j = 0; j < n_trees; j++) {

        const List& phylo_info(chrom_phylo_info[j]);
        std::vector<double> branch_lens = as<std::vector<double>>(
            phylo_info["branch_lens"]);
        arma::Mat<uint64> edges = as<arma::Mat<uint64>>(phylo_info["edges"]);
        std::vector<std::string> tip_labels = as<std::vector<std::string>>(
            phylo_info["labels"]);
        if (branch_lens.size() != edges.n_rows) {
            err_msg = "\nBranch lengths and edges don't have the same ";
            err_msg += "size on chromosome " + std::to_string(i+1) + " and tree ";
            err_msg += std::to_string(j+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }
        if (branch_lens.size() == 0) {
            err_msg = "\nEmpty tree on chromosome " + std::to_string(i+1);
            err_msg += " and tree " + std::to_string(j+1);
            throw(Rcpp::exception(err_msg.c_str(), false));
        }

        n_bases_[j] = as<uint64>(phylo_info["n_bases"]);
        branch_lens_[j] = branch_lens;
        edges_[j] = edges;
        tip_labels_[j] = tip_labels;
    }


    *this = PhyloOneChrom(n_bases_, branch_lens_, edges_, tip_labels_, mutator_base);

    return;

}





/*
 Evolve all trees.
 */
int PhyloOneChrom::evolve(pcg64& eng,
                          Progress& prog_bar) {

    for (uint64 i = 0; i < trees.size(); i++) {

        if (i > 0) {
            for (uint64 j = 0; j < trees[i].ends.size(); j++) {
                trees[i].starts[j] = trees[i-1].ends[j];
                trees[i].ends[j] = trees[i-1].ends[j] + trees[i].end - trees[i-1].end;
            }
        }

#ifdef __JACKALOPE_DIAGNOSTICS
        Rcout << "-- tree " << i << std::endl;
#endif

        int status = one_tree(i, eng, prog_bar);
        if (status < 0) return status;
    }

    return 0;
}






PhyloInfo::PhyloInfo(const List& genome_phylo_info, const TreeMutator& mutator_base) {

    uint64 n_chroms = genome_phylo_info.size();

    if (n_chroms == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }

    phylo_one_chroms = std::vector<PhyloOneChrom>(n_chroms);

    // Fill tree and mutator info (i.e., everything but haplotype info):
    for (uint64 i = 0; i < n_chroms; i++) {
        phylo_one_chroms[i].fill_tree_mutator(genome_phylo_info, i, mutator_base);
    }
}







/*
 Evolve all chromosomes along trees.
*/
XPtr<HapSet> PhyloInfo::evolve_chroms(
        SEXP& ref_genome_ptr,
        const uint64& n_threads,
        const bool& show_progress) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    // (I'm simply extracting tip labels from the first tree, as they should all be
    // the same due to the process_phy function in R/create_haplotypes.R)
    XPtr<HapSet> hap_set(new HapSet(*ref_genome, phylo_one_chroms[0].trees[0].tip_labels),
                         true);

    uint64 n_chroms = ref_genome->size();
    uint64 total_chrom = ref_genome->total_size;

    Progress prog_bar(total_chrom, show_progress);
    std::vector<int> status_codes(n_threads, 0);


    if (n_chroms != phylo_one_chroms.size()) {
        std::string err_msg = "\n# items in phylo. info must be of same length as ";
        err_msg += "# chromosomes in reference genome";
        throw(Rcpp::exception(err_msg.c_str(), false));
    }

    // Generate seeds for random number generators (1 RNG per thread)
    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_threads);

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
{
#endif

    std::vector<uint64> active_seeds;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint64 active_thread = omp_get_thread_num();
#else
    uint64 active_thread = 0;
#endif
    int& status_code(status_codes[active_thread]);
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint64 i = 0; i < n_chroms; i++) {

        if (status_code != 0) continue;

#ifdef __JACKALOPE_DIAGNOSTICS
        Rcout << std::endl << ">> chrom " << i << std::endl;
#endif

        PhyloOneChrom& chrom_phylo(phylo_one_chroms[i]);

        // Set values for haplotype info:
        chrom_phylo.set_hap_info(*hap_set, i);

        // Evolve the chromosome using the chrom_phylo object:
        status_code = chrom_phylo.evolve(eng, prog_bar);

    }

#ifdef _OPENMP
}
#endif

    for (const int& status_code : status_codes) {
        if (status_code == -1) {
            prog_bar.cleanup();
            std::string warn_msg = "\nThe user interrupted phylogenetic evolution. ";
            warn_msg += "Note that changes occur in place, so your haplotypes have ";
            warn_msg += "already been partially added.";
            Rcpp::warning(warn_msg.c_str());
            break;
        }
    }

    return hap_set;

}






//' Evolve all chromosomes in a reference genome.
//'
//' @noRd
//'
//[[Rcpp::export]]
SEXP evolve_across_trees(
        SEXP& ref_genome_ptr,
        const List& genome_phylo_info,
        const std::vector<arma::mat>& Q,
        const std::vector<arma::mat>& U,
        const std::vector<arma::mat>& Ui,
        const std::vector<arma::vec>& L,
        const double& invariant,
        const arma::vec& insertion_rates,
        const arma::vec& deletion_rates,
        const double& epsilon,
        const std::vector<double>& pi_tcag,
        uint64 n_threads,
        const bool& show_progress) {


    // Check that # threads isn't too high and change to 1 if not using OpenMP:
    thread_check(n_threads);


    // Now create mutation sampler:
    TreeMutator mutator(Q, U, Ui, L, invariant,
                        insertion_rates, deletion_rates, epsilon, pi_tcag);


    // Create phylogenetic tree object:
    if (genome_phylo_info.size() == 0) {
        throw(Rcpp::exception("\nEmpty list provided for phylogenetic information.",
                              false));
    }
    PhyloInfo phylo_info(genome_phylo_info, mutator);

    /*
     Now that we have tree(s) and mutator info, we can create haplotypes:
     */

    XPtr<HapSet> hap_set = phylo_info.evolve_chroms(ref_genome_ptr,
                                                    n_threads, show_progress);


    return hap_set;
}







