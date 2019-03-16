

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng
#include <string>  // string class
#include <random>  // distributions



#include "gemino_types.hpp"  // uint32
#include "seq_classes_ref.hpp"  // Ref* classes
#include "seq_classes_var.hpp"  // Var* classes
#include "pcg.hpp"  // runif_01
#include "table_sampler.hpp"  // TableSampler
#include "alias_sampler.hpp"  // AliasSampler
#include "util.hpp"  // clear_memory
#include "str_manip.hpp"  // rev_comp
#include "sequencer.hpp"  // SequenceIdentifierInfo
#include "pacbio.hpp"  // PacBio* types





double PacBioQualityError::calc_min_exp() {

    double min_exp_ = 1;

    double total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
        std::pow(prob_subst, min_exp_);
    double left, right;

    if (total < prob_thresh) {
        while (total < prob_thresh) {
            min_exp_ /= 2;
            total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
                std::pow(prob_subst, min_exp_);
        }
        left = min_exp_;
        right = min_exp_ * 2;
    } else {
        while (total > prob_thresh) {
            min_exp_ *= 2;
            total = std::pow(prob_ins, min_exp_) + std::pow(prob_del, min_exp_) +
                std::pow(prob_subst, min_exp_);
        }
        left = min_exp_ / 2;
        right = min_exp_;
    }
    for (uint32 i = 0; i < 15; i++) {
        double m = (left + right) / 2;
        total = std::pow(prob_ins, m) + std::pow(prob_del, m) + std::pow(prob_subst, m);
        if (total == prob_thresh){
            min_exp_ = m;
            break;
        } else if (total > prob_thresh) {
            left = m;
            min_exp_ = (m + right) / 2;
        } else {
            right = m;
            min_exp_ = (left + m) / 2;
        }
    }

    return min_exp_;
}




void PacBioQualityError::modify_probs(pcg64& eng,
                                      const double& passes_left,
                                      const double& passes_right) {


    // Thresholds for truncated normal distribution
    double left_thresh = (min_exp - (std::sqrt(passes_left + sqrt_params[0]) -
                          sqrt_params[1]))  / sigmoid(passes_left);
    double right_thresh = (min_exp - (std::sqrt(passes_right + sqrt_params[0]) -
                           sqrt_params[1])) / sigmoid(passes_right);

    // Sample for noise:
    double incr_quals_l = trunc_norm(left_thresh, eng);
    double incr_quals_r = trunc_norm(right_thresh, eng);

    // Exponents of quality increase
    double exponent_left = incr_quals_l * sigmoid(passes_left) +
        std::sqrt(passes_left + sqrt_params[0]) - sqrt_params[1];
    double exponent_right = incr_quals_r * sigmoid(passes_right) +
        std::sqrt(passes_right + sqrt_params[0]) - sqrt_params[1];
    // SimLoRD prevents from going below 0.6 by chance:
    if (exponent_left < 0.6) exponent_left = 0.6;
    if (exponent_right < 0.6) exponent_right = 0.6;

    // Now calculate new cumulative error probabilities:
    cum_probs_left[0] = std::pow(prob_ins, exponent_left);
    cum_probs_left[1] = std::pow(prob_del, exponent_left) + cum_probs_left[0];
    cum_probs_left[2] = std::pow(prob_subst, exponent_left) + cum_probs_left[1];

    cum_probs_right[0] = std::pow(prob_ins, exponent_right);
    cum_probs_right[1] = std::pow(prob_del, exponent_right) + cum_probs_right[0];
    cum_probs_right[2] = std::pow(prob_subst, exponent_right) + cum_probs_right[1];

    return;

}



