
#include "mutator_indels.h"  // IndelMutator and debugging preprocessor directives

/*
 This defines classes for adding insertions and deletions.
 */


#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <cmath>  // pow
#include <vector>  // vector class
#include <string>  // string class
#include <random>  // poisson_distribution


#include "var_classes.h"  // Var* classes
#include "pcg.h"  // runif_01()



using namespace Rcpp;



/*
 This calculates...
 - tau (the period of time over which to generate indels)
 - rates over the whole chromosome and `tau` time units (`rates_tau`)
 - new branch length after progressing `tau` time units (`b_len`)
 */
void IndelMutator::calc_tau(double& b_len) {

    const double chrom_size(var_chrom->chrom_size);

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




void IndelMutator::add_indels(double b_len,
                              const uint64& begin,
                              const uint64& end,
                              SubMutator& subs,
                              pcg64& eng) {

#ifdef __JACKALOPE_DEBUG
    if (!var_chrom) stop("var_chrom is nullptr in add_indels");
    if (b_len < 0) stop("b_len < 0 in add_indels");
    if (begin >= var_chrom->size()) stop("begin >= var_chrom->size() in add_indels");
    if (end > var_chrom->size()) stop("end > var_chrom->size() in add_indels");
#endif

    while (b_len > 0) {

        calc_tau(b_len);

        std::deque<uint32> non_zeros;
        uint32 total_events = 0;

        for (uint32 i = 0; i < rates_tau.n_elem; i++) {

            distr.param(std::poisson_distribution<uint32>::param_type(rates_tau(i)));

            n_events[i] = distr(eng);
            if (n_events[i] > 0U) {
                non_zeros.push_back(i);
                total_events += n_events[i];
            }

        }

        // Adding indels in random order:
        while (total_events > 0) {

            uint32 nz_i = static_cast<uint32>(runif_01(eng) * non_zeros.size());
            uint32 ch_i = non_zeros[nz_i];

            if (changes(ch_i) > 0) {
                uint64 size = static_cast<uint64>(changes(ch_i));
                uint64 pos = static_cast<uint64>(runif_01(eng) * (end - begin) + begin);
                std::string str(size, 'x');
                insert.sample(str, eng);
                var_chrom->add_insertion(str, pos);
                subs.insertion_adjust(size, pos, eng);
            } else {
                uint64 size = std::min(static_cast<uint64>(std::abs(changes(ch_i))),
                                       end - begin);
                uint64 pos = static_cast<uint64>(runif_01(eng) *
                    (end - begin - size + 1) + begin);
                var_chrom->add_deletion(size, pos);
                subs.deletion_adjust(size, pos);
            }

            n_events[ch_i]--;
            if (n_events[ch_i] == 0) non_zeros.erase(non_zeros.begin() + ch_i);

            total_events--;
        }


    }

    return;
}


