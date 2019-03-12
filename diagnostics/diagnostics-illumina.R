

#'
#' This tests that the Illumina sequencer produces quality profiles similar
#' to inputs
#'


library(gemino)
library(tidyverse)
source(".Rprofile")

# Profile I'll be using
prof <- gemino:::read_profile(NULL, "HS25", read_length = 100, 1)

# Predicted qualities for each nucleotide
pred_df <- map_dfr(1:4,
                   function(nti) {
                       map_dfr(1:100,
                               function(i) {
                                   cbind(nti, i, prof$quals[[nti]][[i]],
                                         prof$qual_probs[[nti]][[i]]) %>%
                                       as.data.frame() %>%
                                       as_tibble() %>%
                                       set_names(c("nucleo", "pos", "qual", "prob")) %>%
                                       mutate_at(vars(pos, qual), as.integer) %>%
                                       mutate(nucleo = c("T", "C", "A", "G")[nti])
                               })
                   })




# Extracting integer qualities from a string:
Rcpp::sourceCpp(code = "
#include <Rcpp.h>
#include <string>
#include <vector>

//[[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//[[Rcpp::export]]
std::vector<std::vector<int>> quals_by_pos(const std::vector<std::string>& input) {

    size_t n_reads = input.size();
    size_t n_pos = input[0].size();
    std::vector<std::vector<int>> output(n_pos);

    for (size_t i = 0; i < n_pos; i++) {
        output[i].reserve(n_reads);
        for (size_t j = 0; j < n_reads; j++) {
            output[i].push_back(static_cast<int>(input[j][i]) - static_cast<int>('!'));
        }
    }

    return output;
}

/*
 Indicates which nucleotide mostly comprises each read.
 Useful for distinguishing reverse complements.
 */
//[[Rcpp::export]]
std::vector<int> which_nt(const std::vector<std::string>& reads) {

    size_t n_reads = reads.size();
    size_t n_pos = reads[0].size();
    std::vector<int> output(n_reads);
    std::vector<int> nt_map(256, 0);
    std::string nts = \"TCAG\";
    for (int i = 0; i < nts.size(); i++) nt_map[nts[i]] = i;

    std::vector<int> nt_counts(4);

    for (size_t i = 0; i < n_reads; i++) {
        for (int& c : nt_counts) c = 0;
        for (size_t j = 0; j < n_pos; j++) {
            nt_counts[nt_map[reads[i][j]]]++;
        }
        int max_nt_ind = std::max_element(nt_counts.begin(), nt_counts.end()) -
            nt_counts.begin();
        output[i] = max_nt_ind;
    }

    return output;
}

/*
 Mismatches by position (assuming indel probs are zero and
 that reads are mostly the same nucleotide!)
*/
//[[Rcpp::export]]
std::vector<int> mm_by_pos(const std::vector<std::string>& reads) {

    size_t n_reads = reads.size();
    size_t n_pos = reads[0].size();
    std::vector<int> output(n_pos, 0);
    std::vector<int> nt_map(256, 0);
    std::string nts = \"TCAG\";
    for (int i = 0; i < nts.size(); i++) nt_map[nts[i]] = i;

    std::vector<int> nt_counts(4);
    std::string max_nts(n_reads,'A');

    for (size_t i = 0; i < n_reads; i++) {
        for (int& c : nt_counts) c = 0;
        for (size_t j = 0; j < n_pos; j++) {
            nt_counts[nt_map[reads[i][j]]]++;
        }
        int max_nt_ind = std::max_element(nt_counts.begin(), nt_counts.end()) -
            nt_counts.begin();
        max_nts[i] = nts[max_nt_ind];
    }

    for (size_t i = 0; i < n_reads; i++) {
        for (size_t j = 0; j < n_pos; j++) {
            if (reads[i][j] != max_nts[i]) output[j]++;
        }
    }

    return output;
}

//[[Rcpp::export]]
std::vector<int> mm_by_qual(const std::vector<std::string>& reads,
                            const std::vector<std::string>& quals,
                            const int& max_qual) {

    size_t n_reads = reads.size();
    size_t n_pos = reads[0].size();
    std::vector<int> output(max_qual, 0);
    std::vector<int> nt_map(256, 0);
    std::string nts = \"TCAG\";
    for (int i = 0; i < nts.size(); i++) nt_map[nts[i]] = i;

    std::vector<int> nt_counts(4);
    std::string max_nts(n_reads,'A');

    for (size_t i = 0; i < n_reads; i++) {
        for (int& c : nt_counts) c = 0;
        for (size_t j = 0; j < n_pos; j++) {
            nt_counts[nt_map[reads[i][j]]]++;
        }
        int max_nt_ind = std::max_element(nt_counts.begin(), nt_counts.end()) -
            nt_counts.begin();
        max_nts[i] = nts[max_nt_ind];
    }

    for (size_t i = 0; i < n_reads; i++) {
        for (size_t j = 0; j < n_pos; j++) {
            if (reads[i][j] != max_nts[i]) {
                int qual = static_cast<int>(quals[i][j]) - static_cast<int>('!');
                output[qual]++;
            }
        }
    }

    return output;
}
")






# ========================================================
# ========================================================

# Testing mismatches and quality scores:

# ========================================================
# ========================================================


# test_mm_quals <- function(nt) {

bases <- c("T", "C", "A", "G")

# 1 sequence of length 100e3
seq <- sapply(bases, function(nt) paste(rep(nt, 100e3), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

illumina(rg, out_prefix = "~/Desktop/fq/test",
         n_reads = 100e3, read_length = 100, paired = FALSE, seq_sys = "HS25",
         frag_mean = 400, frag_sd = 100, ins_prob1 = 0, del_prob1 = 0)

fq <- readLines("~/Desktop/fq/test_R1.fq")
reads <- fq[seq(2, length(fq), 4)]
quals <- fq[seq(4, length(fq), 4)]
nts <- which_nt(reads) + 1

# Should be ~0.25 each:
map_dbl(1:4, ~ mean(nts == .x))

#'
#' Looking at mapping qualities as predicted by position on read and
#' starting nucleotide
#'

qbp_df <- map_dfr(1:4, function(nt) {
    do.call(cbind, quals_by_pos(quals[nts == nt])) %>%
        as.data.frame() %>%
        as_tibble() %>%
        set_names(sprintf("p%i", 1:100)) %>%
        gather("position", "quality") %>%
        mutate(position = as.integer(gsub("p", "", position))) %>%
        group_by(position, quality) %>%
        summarize(count = n(),
                  pred_count = sum(nts == nt) *
                      pred_df$prob[pred_df$pos == position[1] &
                                       pred_df$qual == quality[1] &
                                       pred_df$nucleo == bases[nt]]) %>%
        ungroup() %>%
        mutate(nucleo = bases[nt])
}) %>%
    mutate(nucleo = factor(nucleo, levels = bases))


lims <- log(range(c(qbp_df$count, qbp_df$pred_count)))
qbp_df %>%
    ggplot(aes(log(pred_count), log(count))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 12)) +
    facet_wrap(~ nucleo, nrow = 2) +
    scale_y_continuous("Observed quality count", breaks = log(40^(0:3)),
                       labels = 40^(0:3), limits = lims) +
    scale_x_continuous("Predicted quality count", breaks = log(40^(0:3)),
                       labels = 40^(0:3), limits = lims)



#'
#' Now looking at mismatches as predicted by quality:
#'

qbp_all <- do.call(c, quals_by_pos(quals))
mmbq_df <- tibble(mm = mm_by_qual(reads, quals, max_qual = 50) %>%
                      keep(~ . > 0),
                  qual = sort(unique(qbp_all)))%>%
    mutate(n = map_dbl(qual, ~ sum(qbp_all == .x)),
           p = 10^(qual / -10.0),
           mm_pred = n * p)


lims <- log(range(c(mmbq_df$mm, mmbq_df$mm_pred)))

mmbq_df %>%
    ggplot(aes(log(mm_pred), log(mm))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted mismatch count", breaks = log(10^(1:3)),
                       labels = 10^(1:3), limits = lims) +
    scale_y_continuous("Observed mismatch count", breaks = log(10^(1:3)),
                       labels = 10^(1:3), limits = lims) +
    theme_classic()




# ========================================================
# ========================================================

# Testing indel rates

# ========================================================
# ========================================================


# Making fake profile for no mismatches
profile_df <- crossing(nucleo = bases, pos = 0:99, qual = c(255L, 1000L)) %>%
    arrange(nucleo, pos, qual)

write_tsv(profile_df, path = "~/Desktop/test_prof.txt", col_names = FALSE)


# 10 sequences of length 100e3
seqs <- replicate(10, paste(rep("TC", 100e3 / 2), collapse = ""))


# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

illumina(rg, out_prefix = "~/Desktop/fq/test",
         n_reads = 100e3, read_length = 100, paired = FALSE,
         profile1 = "~/Desktop/test_prof.txt",
         frag_mean = 400, frag_sd = 100, ins_prob1 = 0.01, del_prob1 = 0.02)


fq <- readLines("~/Desktop/fq/test_R1.fq")
reads <- fq[seq(2, length(fq), 4)]


Rcpp::sourceCpp(code = "
#include <Rcpp.h>
#include <string>
#include <vector>

//[[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//[[Rcpp::export]]
std::vector<int> count_indels(const std::string& read,
                              const char& motif1,
                              const char& motif2) {

    size_t n_reads = input.size();
    size_t n_pos = read.size();
    std::vector<int> output(2, 0);
    auto iter1 = read.begin();
    auto iter2 = read.begin() + 1;

    if (*iter1 != motif1 || *iter2 != motif2) {
        ;
    }

    while (iter2 < (read.end() - 1)) {

        if (*iter1 != motif1 || *iter2 != motif2) {
            ;
        }

        iter1 += 2;
        iter2 += 2;
    }

    return output;
}
")
