
/*
 ********************************************************

 Methods for evolving chromosomes along phylogenies / gene trees

 ********************************************************
 */



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
#include "var_classes.h"  // Var* classes
#include "mutator.h"  // TreeMutator
#include "pcg.h" // pcg sampler types
#include "phylogenomics.h"
#include "util.h"  // thread_check


using namespace Rcpp;







/*
 Process one phylogenetic tree for a single chromosome with no recombination.

 Note that this function should be changed if any of these VarChroms differ from
 each other (within the range specified if recombination = true).
 They can already have mutations, but to start out, they must all be the same.
 */
int PhyloOneChrom::one_tree(PhyloTree& tree,
                            pcg64& eng,
                            Progress& prog_bar) {

    // For when/if user interrupts function:
    int status;

    // Reset tree of samplers and VarChrom objects representing nodes and tips:
    status = reset(tree, eng, prog_bar);
    if (status < 0) return status;


    uint64 b1, b2;
    VarChrom* chrom1;
    VarChrom* chrom2;
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

        // VarChrom object that parent node refers to:
        chrom1 = &node_chroms[(b1 - n_tips)];
        // For whether chrom we're changing is a tip:
        bool at_tip = b2 < n_tips;
        // Pointer to chrom we're changing:
        chrom2 = (at_tip) ? tip_chroms[b2] : &node_chroms[(b2 - n_tips)];

        // Update rate indices:
        rates[b2] = rates[b1];

        /*
         Update VarChrom objects for this branch.
         Because node VarChrom objects are emptied of mutations for each new tree,
         the below works to add just the mutations related to this tree:
         */
        (*chrom2) += (*chrom1);

        /*
         Now mutate along branch length:
         */
        b_len = tree.branch_lens[i];
        status = mutator.mutate(b_len, *chrom2, eng, prog_bar,
                                tree.start, tree.ends[b2], rates[b2]);
        if (status < 0) return status;


        /*
         To free up some memory, clear info from VarChrom object at `b1` if it's no
         longer needed.
         */
        clear_branches(b1, i, tree);

    }


    // Update progress bar:
    prog_bar.increment(tree.end - tree.start + 1);

    return 0;

}





/*
 Evolve all chromosomes along trees.
*/
XPtr<VarSet> PhyloInfo::evolve_chroms(
        SEXP& ref_genome_ptr,
        const uint64& n_threads,
        const bool& show_progress) {

    XPtr<RefGenome> ref_genome(ref_genome_ptr);

    // (I'm simply extracting tip labels from the first tree, as they should all be
    // the same due to the process_phy function in R/create_variants.R)
    XPtr<VarSet> var_set(new VarSet(*ref_genome, phylo_one_chroms[0].trees[0].tip_labels),
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

        PhyloOneChrom& chrom_phylo(phylo_one_chroms[i]);

        // Set values for variant info:
        chrom_phylo.set_var_info(*var_set, i);

        // Evolve the chromosome using the chrom_phylo object:
        status_code = chrom_phylo.evolve(eng, prog_bar);

    }

#ifdef _OPENMP
}
#endif

    for (const int& status_code : status_codes) {
        if (status_code == -1) {
            std::string warn_msg = "\nThe user interrupted phylogenetic evolution. ";
            warn_msg += "Note that changes occur in place, so your variants have ";
            warn_msg += "already been partially added.";
            Rcpp::warning(warn_msg.c_str());
            break;
        }
    }

    return var_set;

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
     Now that we have tree(s) and mutator info, we can create variants:
     */

    XPtr<VarSet> var_set = phylo_info.evolve_chroms(ref_genome_ptr,
                                                    n_threads, show_progress);


    return var_set;
}







