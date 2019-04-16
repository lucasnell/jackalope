#' jackalope: An efficient, flexible molecular evolution and sequencing simulator.
#'
#' `jackalope` efficiently (i) reads and simulates reference genomes;
#' (ii) generates variants using summary statistics, phylogenies, Variant
#' Call Format (VCF) files, and coalescent simulations—the latter of which can include
#' selection, recombination, and demographic fluctuations;
#' (iii) simulates sequencing error, mapping qualities, and optical/PCR duplicates; and
#' (iv) writes outputs to standard file formats.
#' `jackalope` can simulate single, paired-end, or mate-pair Illumina reads, as well as
#' reads from Pacific BioSciences.
#'
#'
#' @importFrom Rcpp evalCpp
#' @import Rhtslib
#' @useDynLib jackalope, .registration = TRUE
#'
#' @docType package
#' @name jackalope
NULL
