

#'
#' This tests that the PacBio sequencer produces quality profiles similar
#' to inputs
#'


library(gemino)
library(tidyverse)
source(".Rprofile")

Rcpp::sourceCpp("diagnostics/diagnostics-pacbio.cpp")

dir <- paste0(tempdir(check = TRUE), "/")


# 4 sequences of length 100e3
seq <- sapply(c("T", "C", "A", "G"), function(nt) paste(rep(nt, 100e3), collapse = ""))

# Make ref_genome object from a pointer to a RefGenome object based on `seqs`
rg <- ref_genome$new(gemino:::make_ref_genome(seq))

pacbio(rg, out_prefix = paste0(dir, "test"), n_reads = 10e3)

fq <- readLines(paste0(dir, "test_R1.fq"))


# Looking at default read length distribution:

# Defaults from SimLoRD:
test_lens <- function(N) {
    pars <- list(sigma = 0.200110276521, loc = -10075.4363813, scale = 17922.611306,
                 min_len = 50)
    z <- rlnorm(N, log(pars$scale), pars$sigma)
    z <- z + pars$loc
    z[z < pars$min_len] <- pars$min_len
    return(z)
}

read_lens <- nchar(fq[seq(2, length(fq), 4)])

par(mfrow = c(2, 1))
hist(read_lens, main = "Observed read length distribution")
hist(test_lens(10e3), main = "Expected read length distribution")


# Looking at qualities:

quals <- fq[seq(4, length(fq), 4)]
# Mean quality by read:
m_quals <- mean_quals(quals)

# This looks approximately like the SimLoRD panel in Fig. 1 in the paper,
# so I'm pretty satisfied.
# The extra reads of short length here are likely a result of either the way they plotted
# things in the original paper or that they used the genome of *Neurospora crassa*
# rather than simulated reads.


tibble(qual = m_quals, len = read_lens) %>%
    ggplot(aes(len, qual)) +
    geom_hex(bins = 50) +
    xlab("Read length") +
    ylab("Read quality\n(average base quality)") +
    theme_classic() +
    theme(legend.position = c(0.75, 0.5))



# Looking at non-default read length distribution

rl_mat <- cbind(seq(100, 1e4, 100), runif(100))
rl_mat[,2] <- rl_mat[,2] / sum(rl_mat[,2])
pacbio(rg, out_prefix = paste0(dir, "test"), n_reads = 10e3,
       custom_read_lengths = rl_mat)

fq <- readLines(paste0(dir, "test_R1.fq"))
read_lens <- nchar(fq[seq(2, length(fq), 4)])


# Plotting observed vs predicted:

ndrl_df <- tibble(obs = as.integer(table(read_lens)),
                  pred = length(read_lens) * rl_mat[,2])

lims <- range(c(ndrl_df$obs, ndrl_df$pred))

ndrl_df %>%
    ggplot(aes((pred), (obs))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "firebrick2") +
    scale_x_continuous("Predicted read length count",
                       limits = lims) +
    scale_y_continuous("Observed read length count",
                       limits = lims) +
    theme_classic()


