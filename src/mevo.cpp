
#include <RcppArmadillo.h>
#include <cmath>  // pow, log, exp
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>  // vector class
#include <string>  // string class
#include <progress.hpp>  // for the progress bar



#include "gemino_types.h"
#include "mevo.h"
#include "sequence_classes.h"  // Var* and Ref* classes
#include "pcg.h"  // pcg seeding
#include "table_sampler.h"  // table method of sampling
#include "weighted_reservoir.h"  // weighted reservoir sampling
#include "mevo_gammas.h"  // SequenceGammas class
#include "mevo_rate_matrices.h"  // rate matrix functions

using namespace Rcpp;





//' Fill in vectors of mutation probabilities and lengths.
//'
//' These vectors should be initialized already, but there's no need to resize them.
//'
//'
//' @param Q A matrix of substitution rates for each nucleotide.
//' @param xi Overall rate of indels.
//' @param psi Proportion of insertions to deletions.
//' @param pi_tcag Vector of nucleotide equilibrium frequencies for
//'     "T", "C", "A", and "G", respectively.
//' @param rel_insertion_rates Relative insertion rates.
//' @param rel_deletion_rates Relative deletion rates.
//'
//' @noRd
//'
void fill_mut_prob_length_vectors(
        std::vector<std::vector<double>>& probs,
        std::vector<sint32>& mut_lengths,
        const arma::mat& Q,
        const double& xi,
        const double& psi,
        const std::vector<double>& pi_tcag,
        arma::vec rel_insertion_rates,
        arma::vec rel_deletion_rates) {

    // If overall rate of indels is zero, remove all elements from these vectors:
    if (xi <= 0) {
        rel_insertion_rates.reset();
        rel_deletion_rates.reset();
    }
    uint32 n_ins = rel_insertion_rates.n_elem;
    uint32 n_del = rel_deletion_rates.n_elem;
    uint32 n_muts = 4 + n_ins + n_del;

    if (n_muts == 4 && xi > 0) {
        stop("If indel rate > 0, vectors of the relative rates of insertions and "
                 "deletions cannot both be of length 0.");
    }

    // 1 vector of probabilities for each nucleotide: T, C, A, then G
    probs.resize(4);

    if (n_ins > 0) {
        // make relative rates sum to 1:
        rel_insertion_rates /= arma::accu(rel_insertion_rates);
        // Now make them sum to the overall insertion rate:
        double xi_i = (0.25 * xi) / (1 + 1/psi);  // overall insertion rate
        rel_insertion_rates *= xi_i;
    }
    // Same for deletions
    if (n_del > 0) {
        rel_deletion_rates /= arma::accu(rel_deletion_rates);
        double xi_d = (0.25 * xi) / (1 + psi);    // overall deletion rate
        rel_deletion_rates *= xi_d;
    }


    /*
     (1) Combine substitution, insertion, and deletion rates into a single vector
     (2) Create TableSampler for each nucleotide
     (3) Fill the `rates` field with mutation rates for each nucleotide
     */
    for (uint32 i = 0; i < 4; i++) {

        std::vector<double>& qc(probs[i]);

        qc = arma::conv_to<std::vector<double>>::from(Q.row(i));
        // Get the overall mutation rate for this nucleotide
        double qi = -1 * qc[i];
        // Add insertions, then deletions
        qc.reserve(n_muts);
        for (uint32 j = 0; j < rel_insertion_rates.n_elem; j++) {
            qc.push_back(rel_insertion_rates(j));
        }
        for (uint32 j = 0; j < rel_deletion_rates.n_elem; j++) {
            qc.push_back(rel_deletion_rates(j));
        }
        // Divide all by qi to make them probabilities
        for (uint32 j = 0; j < n_muts; j++) qc[j] /= qi;

        // Change the diagonal back to the mutation rate, which will be used later.
        qc[i] = qi;
    }

    // Now filling in mut_lengths vector
    mut_lengths = std::vector<sint32>(n_muts, 0);
    for (uint32 i = 0; i < rel_insertion_rates.n_elem; i++) {
        mut_lengths[i + 4] = static_cast<sint32>(i+1);
    }
    for (uint32 i = 0; i < rel_deletion_rates.n_elem; i++) {
        sint32 ds = static_cast<sint32>(i + 1);
        ds *= -1;
        mut_lengths[i + 4 + rel_insertion_rates.n_elem] = ds;
    }

    return;
}



MutationSampler make_mutation_sampler(VarSequence& var_seq,
                                      const std::vector<std::vector<double>>& probs,
                                      const std::vector<sint32>& mut_lengths,
                                      const std::vector<double>& pi_tcag,
                                      const arma::mat& gamma_mat) {

    MutationTypeSampler mts(probs, mut_lengths);
    TableStringSampler<std::string> tss(mevo::bases, pi_tcag);

    SequenceGammas gammas(gamma_mat);
    std::vector<double> q_tcag(4);
    for (uint32 i = 0; i < 4; i++) q_tcag[i] = probs[i][i];
    MutationRates mr(var_seq, q_tcag, gammas);
    LocationSampler ls(mr);

    MutationSampler ms(var_seq, ls, mts, tss);

    return ms;
}


ChunkMutationSampler make_mutation_sampler(VarSequence& var_seq,
                                           const std::vector<std::vector<double>>& probs,
                                           const std::vector<sint32>& mut_lengths,
                                           const std::vector<double>& pi_tcag,
                                           const arma::mat& gamma_mat,
                                           const uint32& chunk_size) {

    MutationTypeSampler mts(probs, mut_lengths);
    TableStringSampler<std::string> tss(mevo::bases, pi_tcag);

    SequenceGammas gammas(gamma_mat);
    std::vector<double> q_tcag(4);
    for (uint32 i = 0; i < 4; i++) q_tcag[i] = probs[i][i];
    MutationRates mr(var_seq, q_tcag, gammas);
    ChunkLocationSampler ls(mr, chunk_size);

    ChunkMutationSampler ms(var_seq, ls, mts, tss);

    return ms;
}

//' Creates MutationSampler without any of the pointers.
//'
//' `T` should be MutationSampler or ChunkMutationSampler
//' `T` should be LocationSampler or ChunkLocationSampler
//' MutationSampler should always go with LocationSampler, and
//' ChunkMutationSampler with ChunkLocationSampler
//'
//' Before actually using the object output from this function, make sure to...
//' * use `[Chunk]MutationSampler.fill_ptrs(VarSequence& var_seq)` to fill pointers.
//' * use `[Chunk]MutationSampler.fill_gamma(const arma::mat& gamma_mat)` to fill
//'   the gamma matrix.
//' * use `ChunkMutationSampler.location.change_chunk(chunk_size)` if using chunked
//'   version.
//'
//' @noRd
//'
template <typename T, typename U>
XPtr<T> make_mutation_sampler_base_(const arma::mat& Q,
                                    const double& xi,
                                    const double& psi,
                                    const std::vector<double>& pi_tcag,
                                    const arma::vec& rel_insertion_rates,
                                    const arma::vec& rel_deletion_rates) {

    std::vector<std::vector<double>> probs;
    std::vector<sint32> mut_lengths;

    fill_mut_prob_length_vectors(probs, mut_lengths, Q, xi, psi, pi_tcag,
                                 rel_insertion_rates, rel_deletion_rates);

    XPtr<T> out(new T());

    out->type = MutationTypeSampler(probs, mut_lengths);
    out->insert = TableStringSampler<std::string>(mevo::bases, pi_tcag);

    std::vector<double> q_tcag(4);
    for (uint32 i = 0; i < 4; i++) q_tcag[i] = probs[i][i];
    MutationRates mr(q_tcag);
    out->location = U(mr);

    return out;
}

// Wrapper to make non-chunked version available from R

//[[Rcpp::export]]
SEXP make_mutation_sampler_base(const arma::mat& Q,
                                const double& xi,
                                const double& psi,
                                const std::vector<double>& pi_tcag,
                                const arma::vec& rel_insertion_rates,
                                const arma::vec& rel_deletion_rates) {

    XPtr<MutationSampler> out =
        make_mutation_sampler_base_<MutationSampler,LocationSampler>(
                Q, xi, psi, pi_tcag, rel_insertion_rates, rel_deletion_rates);

    return out;
}

// Same thing, but with chunks

//[[Rcpp::export]]
SEXP make_mutation_sampler_chunk_base(const arma::mat& Q,
                                      const double& xi,
                                      const double& psi,
                                      const std::vector<double>& pi_tcag,
                                      const arma::vec& rel_insertion_rates,
                                      const arma::vec& rel_deletion_rates,
                                      const uint32& chunk_size) {

    XPtr<ChunkMutationSampler> out =
        make_mutation_sampler_base_<ChunkMutationSampler,ChunkLocationSampler>(
                Q, xi, psi, pi_tcag, rel_insertion_rates, rel_deletion_rates);

    out->location.change_chunk(chunk_size);

    return out;
}
