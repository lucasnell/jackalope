#' gemino: An efficient, flexible molecular evolution and sequencing simulator.
#'
#' `gemino` efficiently (i) reads and simulates reference genomes;
#' (ii) generates variants using summary statistics, phylogenies, Variant
#' Call Format (VCF) files, and coalescent simulations—the latter of which can include
#' selection, recombination, and demographic fluctuations;
#' (iii) simulates sequencing error, mapping qualities, restriction-enzyme digestion,
#' and variance in coverage among sites; and
#' (iv) writes outputs to standard file formats.
#' `gemino` can simulate single or paired-ended reads for WGS on the Illumina platform,
#' and can be extended to simulate different methods (e.g., original RADseq,
#' double-digest RADseq, and genotyping-by-sequencing), and sequencing technologies
#' (e.g., Pacific BioSciences, Oxford Nanopore Technologies).
#'
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib gemino, .registration = TRUE
#'
#' @docType package
#' @name gemino
NULL
