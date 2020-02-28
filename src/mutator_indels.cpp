

/*
 This defines classes for adding insertions and deletions.
 */


#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <cmath>  // pow
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // poisson_distribution


#include "mutator_indels.h"  // IndelMutator
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // runif_01()
#include "util.h"  // interrupt_check


using namespace Rcpp;



// Add just one indel
inline void IndelMutator::one_indel__(const double& change,
                                      std::string& insert_str,
                                      const uint64& begin,
                                      uint64& end,
                                      std::deque<uint8>& rate_inds,
                                      SubMutator& subs,
                                      HapChrom& hap_chrom,
                                      pcg64& eng) {

    if (change > 0) {
        uint64 size = static_cast<uint64>(change);
        uint64 pos = static_cast<uint64>(runif_01(eng) * (end - begin) + begin);
        insert_str.clear();
        for (uint32 j = 0; j < size; j++) insert_str += insert.sample(eng);
        hap_chrom.add_insertion(insert_str, pos);
        subs.insertion_adjust(size, pos, begin, rate_inds, eng);
        end += size;
    } else {
        uint64 size = std::min(static_cast<uint64>(std::abs(change)),
                               end - begin);
        uint64 pos = static_cast<uint64>(runif_01(eng) *
            (end - begin - size + 1) + begin);
        hap_chrom.add_deletion(size, pos);
        subs.deletion_adjust(size, pos, begin, rate_inds);
        end -= size;
    }

    return;

}



inline void IndelMutator::exact_sim(int& status,
                                    double& b_len,
                                    const uint64& begin,
                                    uint64& end,
                                    std::deque<uint8>& rate_inds,
                                    SubMutator& subs,
                                    HapChrom& hap_chrom,
                                    pcg64& eng,
                                    Progress& prog_bar) {


    // For insertions:
    std::string insert_str;
    insert_str.reserve(rates.n_elem / 2);

    uint32 iters = 0;

    double rate = total_rate * (end - begin);

    jump_distr.param(std::exponential_distribution<double>::param_type(rate));

    b_len -= jump_distr(eng);

    uint32 i; // index for which indel-event type


#ifdef __JACKALOPE_DIAGNOSTICS
    double csize;
    arma::vec n_muts;
#endif

    while (b_len > 0) {

        if (interrupt_check(iters, prog_bar)) {
            status = -1;
            return;
        }

        i = event_sampler.sample(eng);

        one_indel__(changes(i), insert_str, begin, end, rate_inds, subs, hap_chrom, eng);

        if (end == begin) return;

        rate = total_rate * (end - begin);
        jump_distr.param(std::exponential_distribution<double>::param_type(rate));
        b_len -= jump_distr(eng);

    }


    return;
}

/*
 This calculates...
 - tau (the period of time over which to generate indels)
 - rates over the whole chromosome and `tau` time units (`rates_tau`)
 - new branch length after progressing `tau` time units (`b_len`)
 */
void IndelMutator::calc_tau(double& b_len, HapChrom& hap_chrom) {

    const double chrom_size(hap_chrom.chrom_size);

    // Now rates are in units of "indels per unit time" (NOT yet over `tau` time units):
    rates_tau = rates * chrom_size;

    // For the expected number of bp changes per time over entire chromosome...
    double mu = arma::accu(changes % rates_tau);             // mean
    double sig = arma::accu(changes % changes % rates_tau);  // variance

    tau = std::min(std::max(eps * chrom_size, 1.0) / std::abs(mu),
                   std::pow(std::max(eps * chrom_size, 1.0), 2U) / sig);

    // We don't want to exceed the remaining branch length
    if (b_len < tau) tau = b_len;

    // Adjust the remaining branch length
    b_len -= tau;

    // Now this vector is in units of "indels per `tau` time units"
    rates_tau *= tau;

    return;

}



// Add indels, adjust `end` (`end == begin` when chromosome region is of size zero)
int IndelMutator::add_indels(double b_len,
                             const uint64& begin,
                             uint64& end,
                             std::deque<uint8>& rate_inds,
                             SubMutator& subs,
                             HapChrom& hap_chrom,
                             pcg64& eng,
                             Progress& prog_bar) {

#ifdef __JACKALOPE_DEBUG
    if (b_len < 0) {
        Rcout << std::endl << b_len << std::endl;
        stop("b_len < 0 in add_indels");
    }
    if (begin >= hap_chrom.size()) {
        Rcout << std::endl << begin << ' ' << hap_chrom.size() << std::endl;
        stop("begin >= hap_chrom.size() in add_indels");
    }
    if (end > hap_chrom.size()) {
        Rcout << std::endl << end << ' ' << hap_chrom.size() << std::endl;
        stop("end > hap_chrom.size() in add_indels");
    }
#endif

    if ((rates.n_elem == 0) || (b_len == 0) || (end == begin)) return 0;
    if (prog_bar.is_aborted() || prog_bar.check_abort()) return -1;

    // Doing exact simulation in this case
    if (eps <= 0) {
        int status = 0;
        exact_sim(status, b_len, begin, end, rate_inds, subs, hap_chrom, eng, prog_bar);
        return status;
    }

    // Vector of indel-type indices, one item per indel "event"
    std::vector<uint32> events;
    // For insertions:
    std::string insert_str;
    insert_str.reserve(rates.n_elem / 2);

    uint32 iters = 0;

#ifdef __JACKALOPE_DIAGNOSTICS
    double csize;
    arma::vec n_muts;
#endif

    while (b_len > 0) {

        /*
         ----------------
         Determine how many of each indel-type occur over `tau` time units:
         ----------------
         */

        calc_tau(b_len, hap_chrom);

#ifdef __JACKALOPE_DIAGNOSTICS
        csize = hap_chrom.size();
        n_muts.zeros(rates_tau.n_elem);
#endif

        // Reset `events` between rounds:
        if (events.size() > 0) events.clear();
        // Reserve memory (`arma::accu(rates_tau)` is the expected value of the # events)
        events.reserve(1 + 1.5 * arma::accu(rates_tau));

        for (uint32 i = 0; i < rates_tau.n_elem; i++) {

            distr.param(std::poisson_distribution<uint32>::param_type(rates_tau(i)));

            uint32 n_events = distr(eng);

#ifdef __JACKALOPE_DIAGNOSTICS
            n_muts(i) += n_events;
#endif

            for (uint32 j = 0; j < n_events; j++) events.push_back(i);

            if (interrupt_check(iters, prog_bar, 10)) return -1;

        }

        /*
         ----------------
         Adding indels in random order:
         ----------------
         */
        jlp_shuffle<std::vector<uint32>>(events, eng);   // shuffle them first

        iters = 0;

        for (uint32 i = 0; i < events.size(); i++) {

            // The amount that this indel-type changes the chromosome size:
            double& change(changes(events[i]));

#ifdef __JACKALOPE_DEBUG
            if (change == 0) stop("change == 0 inside add_indels");
#endif

            // Check for user interrupt every 1000 indels:
            if (interrupt_check(iters, prog_bar)) return -1;

            one_indel__(change, insert_str, begin, end, rate_inds, subs, hap_chrom, eng);

            if (end == begin) return 0;

        }

#ifdef __JACKALOPE_DIAGNOSTICS
        Rcout << "+- " << csize << ' ' << tau << " | ";
        for (double& n : n_muts) Rcout << n << ' ';
        Rcout << std::endl;
#endif

    }


    return 0;
}


