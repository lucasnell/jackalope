#ifndef __JACKALOPE_MUTATOR_SUBS_H
#define __JACKALOPE_MUTATOR_SUBS_H



/*
 This defines classes for adding substitutions.
 */

#include "jackalope_config.h" // controls debugging and diagnostics output

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <progress.hpp>  // for the progress bar
#include <vector>  // vector class
#include <string>  // string class


#include "jackalope_types.h" // integer types
#include "hap_classes.h"  // Hap* classes
#include "pcg.h"  // pcg seeding
#include "alias_sampler.h"  // alias method of sampling
#include "util.h"  // str_stop



using namespace Rcpp;


// All return 4 except for TCAG
inline std::vector<uint8> make_char_map() {
    std::vector<uint8> out(256, 4);
    std::string bases = "TCAG";
    for (uint32 i = 0; i < 4; i++) out[bases[i]] = i;
    return out;
}







class SubMutator {

public:

    std::vector<arma::mat> Q;
    std::vector<arma::mat> U;
    std::vector<arma::mat> Ui;
    std::vector<arma::vec> L;
    double invariant;
    const std::vector<uint8> char_map = make_char_map();
    std::vector<std::vector<AliasSampler>> samplers;
    std::vector<arma::mat> Pt;


    SubMutator() {}
    SubMutator(const std::vector<arma::mat>& Q_,
               const std::vector<arma::mat>& U_,
               const std::vector<arma::mat>& Ui_,
               const std::vector<arma::vec>& L_,
               const double& invariant_)
        : Q(Q_), U(U_), Ui(Ui_), L(L_), invariant(invariant_),
          samplers(Q_.size(), std::vector<AliasSampler>(4)),
          Pt(Q_.size(), arma::mat(4,4)),
          site_var(((invariant_ > 0) || (Q_.size() > 1)) ? true : false) {
#ifdef __JACKALOPE_DEBUG
        if (Q_.size() == 0) stop("in SubMutator constr, Q_.size() == 0");
        if (Q_.size() > 255) stop("in SubMutator constr, Q_.size() > 255");
        if (U_.size() > 255) stop("in SubMutator constr, U_.size() > 255");
        if (Ui_.size() > 255) stop("in SubMutator constr, Ui_.size() > 255");
        if (L_.size() > 255) stop("in SubMutator constr, L_.size() > 255");
#endif
    }


    SubMutator(const SubMutator& other)
        : Q(other.Q), U(other.U), Ui(other.Ui), L(other.L), invariant(other.invariant),
          samplers(other.samplers), Pt(other.Pt),
          site_var(other.site_var) {};

    SubMutator& operator=(const SubMutator& other) {
        Q = other.Q;
        U = other.U;
        Ui = other.Ui;
        L = other.L;
        invariant = other.invariant;
        samplers = other.samplers;
        Pt = other.Pt;
        site_var = other.site_var;
        return *this;
    }


    int new_rates(const uint64& begin,
                  const uint64& end,
                  std::deque<uint8>& rate_inds,
                  pcg64& eng,
                  Progress& prog_bar);

    int add_subs(const double& b_len,
                 const uint64& begin,
                 const uint64& end,
                 const std::deque<uint8>& rate_inds,
                 HapChrom& hap_chrom,
                 pcg64& eng,
                 Progress& prog_bar);

    // Adjust rate_inds for indels:
    void deletion_adjust(const uint64& size, uint64 pos, const uint64& begin,
                         std::deque<uint8>& rate_inds);
    void insertion_adjust(const uint64& size, uint64 pos, const uint64& begin,
                          std::deque<uint8>& rate_inds, pcg64& eng);


private:

    bool site_var; // for whether to include among-site variability

    inline void adjust_mats(const double& b_len);

    inline void subs_before_muts__(const uint64& pos,
                                   uint64& mut_i,
                                   const std::string& bases,
                                   const uint8& rate_i,
                                   HapChrom& hap_chrom,
                                   pcg64& eng);
    inline int subs_before_muts(const uint64& begin,
                                const uint64& end,
                                uint64& mut_i,
                                const uint8& max_gamma,
                                const std::string& bases,
                                const std::deque<uint8>& rate_inds,
                                HapChrom& hap_chrom,
                                pcg64& eng,
                                Progress& prog_bar,
                                uint32& iters);

    inline void subs_after_muts__(const uint64& pos,
                                  uint64& mut_i,
                                  const std::string& bases,
                                  const uint8& rate_i,
                                  HapChrom& hap_chrom,
                                  pcg64& eng);
    inline int subs_after_muts(uint64& pos,
                               const uint64& begin,
                               const uint64& end1,
                               const uint64& end2,
                               uint64& mut_i,
                               const uint8& max_gamma,
                               const std::string& bases,
                               const std::deque<uint8>& rate_inds,
                               HapChrom& hap_chrom,
                               pcg64& eng,
                               Progress& prog_bar,
                               uint32& iters);



};






//' Changing P(t) matrix with new branch lengths or times.
//'
//'
//' Equivalent to `U %*% diag(exp(L * t)) %*% Ui`
//'
//' @noRd
//'
inline void Pt_calc(const arma::mat& U,
                    const arma::mat& Ui,
                    const arma::vec& L,
                    const double& t,
                    arma::mat& Pt) {

    arma::mat diag_L = arma::diagmat(arma::exp(L * t));

    Pt = U * diag_L * Ui;

    return;
}

//' Calculating P(t) using repeated matrix squaring, for UNREST model only.
//'
//' @noRd
//'
inline void Pt_calc(const arma::mat& Q,
                    const uint32& k,
                    const double& t,
                    arma::mat& Pt) {

    double m = static_cast<double>(1U<<k);

    Pt = arma::eye<arma::mat>(4, 4) + Q * t / m + 0.5 * (Q * t / m) * (Q * t / m);

    for (uint32 i = 0; i < k; i++) Pt = Pt * Pt;

    return;

}










#endif
