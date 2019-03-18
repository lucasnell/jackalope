

#'
#' This tests that the Illumina sequencer produces quality profiles similar
#' to inputs
#'


library(gemino)
library(tidyverse)
source(".Rprofile")

Rcpp::sourceCpp("diagnostics/diagnostics-illumina.cpp")

dir <- paste0(tempdir(check = TRUE), "/")

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


# Making fake profile for no mismatches
profile_df <- crossing(nucleo = c("T", "C", "A", "G"),
                       pos = 0:99,
                       qual = c(255L, 1000L)) %>%
    arrange(nucleo, pos, qual)

write_tsv(profile_df, path = paste0(dir, "test_prof.txt"), col_names = FALSE)




# ========================================================
# ========================================================

# Testing mismatches and quality scores:

# ========================================================
# ========================================================


bases <- c("T", "C", "A", "G")

# 1 sequence of length 100e3
seq <- sapply(bases, function(nt) paste(rep(nt, 100e3), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

illumina(rg, out_prefix = paste0(dir, "test"),
         n_reads = 100e3, read_length = 100, paired = FALSE, seq_sys = "HS25",
         frag_mean = 400, frag_sd = 100, ins_prob1 = 0, del_prob1 = 0)

fq <- readLines(paste0(dir, "test_R1.fq"))
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



# 100e3 sequences of length 10
seqs <- replicate(100e3, paste(rep("TA", 10 / 2), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
rg <- ref_genome$new(gemino:::make_ref_genome(seqs))



indel_test <- function(.ins_prob, .del_prob) {
    stopifnot(.ins_prob == 0 | .del_prob == 0)
    stopifnot(.ins_prob > 0 | .del_prob > 0)
    illumina(rg, out_prefix = paste0(dir, "test"),
             n_reads = 100e3, read_length = 100, paired = FALSE,
             profile1 = paste0(dir, "test_prof.txt"),
             frag_mean = 400, frag_sd = 100,
             ins_prob1 = .ins_prob,
             del_prob1 = .del_prob)
    fq <- readLines(paste0(dir, "test_R1.fq"))
    .lens <- nchar(fq[seq(2, length(fq), 4)])
    if (.ins_prob > 0) {
        .type = "insertion"
        .n <- sum(.lens - 10)
        .n_pred <- .ins_prob * 10 * 100e3
    } else if (.del_prob > 0) {
        .type = "deletion"
        .n <- sum(10 - .lens)
        .n_pred <- .del_prob * 10 * 100e3
    }
    return(tibble(type = .type, n = .n, n_pred = .n_pred))
}




# 20 combos takes ~ 5 sec
indel_df <- crossing(ins_prob = c(0, 10^(seq(-1,-4,length.out = 11))),
                     del_prob = c(0, 10^(seq(-1,-4,length.out = 11)))) %>%
    filter((ins_prob == 0 | del_prob == 0) & (ins_prob > 0 | del_prob > 0)) %>%
    pmap_dfr(~ indel_test(..1, ..2))

lims <- log(range(c(indel_df$n, indel_df$n_pred)))
labs <- 10^(2:5)

indel_df %>%
    ggplot(aes(log(n_pred), log(n))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick") +
    geom_point(aes(color = type), alpha = 0.5, size = 3) +
    scale_x_continuous("Predicted indel count", limits = lims,
                       breaks = log(labs), labels = labs) +
    scale_y_continuous("Observed indel count", limits = lims,
                       breaks = log(labs), labels = labs) +
    scale_color_brewer(NULL, palette = "Dark2") +
    theme_classic() +
    theme(legend.position = c(0.1, 0.9)) +
    NULL





# ========================================================
# ========================================================

# Testing paired-end reads:

# ========================================================
# ========================================================


# 1 sequence of length 200
seq <- paste(c(rep('C', 25), rep('N', 150), rep('T', 25)), collapse = "")

poss_pairs <- c(paste(c(rep('C', 25), rep('N', 75)), collapse = ""),
                paste(c(rep('A', 25), rep('N', 75)), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seq`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

illumina(rg, out_prefix = paste0(dir, "test"),
         n_reads = 10e3, read_length = 100,
         # Paired-end reads:
         paired = TRUE, matepair = FALSE,
         # Fragments will always be of length 200:
         frag_mean = 400, frag_sd = 100,
         frag_len_min = 200, frag_len_max = 200,
         # No sequencing errors:
         ins_prob1 = 0, del_prob1 = 0,
         ins_prob2 = 0, del_prob2 = 0,
         profile1 = paste0(dir, "test_prof.txt"),
         profile2 = paste0(dir, "test_prof.txt"))

fq1 <- readLines(paste0(dir, "test_R1.fq"))
fq2 <- readLines(paste0(dir, "test_R2.fq"))
reads1 <- fq1[seq(2, length(fq1), 4)]
reads2 <- fq2[seq(2, length(fq2), 4)]

# Should both be true:
all(unique(reads1) %in% poss_pairs)
all(unique(reads2) %in% poss_pairs)







# ========================================================
# ========================================================

# Testing mate-pair reads:

# ========================================================
# ========================================================

# 1 sequence of length 200
seq <- paste(c(rep('C', 25), rep('N', 150), rep('T', 25)), collapse = "")

poss_pairs <- c(paste(c(rep('N', 75), rep('T', 25)), collapse = ""),
                paste(c(rep('N', 75), rep('G', 25)), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seq`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

illumina(rg, out_prefix = paste0(dir, "test"),
         n_reads = 10, read_length = 100,
         # Mate-pair reads:
         paired = TRUE, matepair = TRUE,
         # Fragments will always be of length 200:
         frag_mean = 400, frag_sd = 100,
         frag_len_min = 200, frag_len_max = 200,
         # No sequencing errors:
         ins_prob1 = 0, del_prob1 = 0,
         ins_prob2 = 0, del_prob2 = 0,
         profile1 = paste0(dir, "test_prof.txt"),
         profile2 = paste0(dir, "test_prof.txt"))

fq1 <- readLines(paste0(dir, "test_R1.fq"))
fq2 <- readLines(paste0(dir, "test_R2.fq"))
reads1 <- fq1[seq(2, length(fq1), 4)]
reads2 <- fq2[seq(2, length(fq2), 4)]

# Should both be true:
all(unique(reads1) %in% poss_pairs)
all(unique(reads2) %in% poss_pairs)

