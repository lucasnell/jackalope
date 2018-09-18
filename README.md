
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build
Status](https://travis-ci.com/lucasnell/gemino.svg?branch=master)](https://travis-ci.com/lucasnell/gemino)
[![codecov](https://codecov.io/gh/lucasnell/gemino/branch/master/graph/badge.svg)](https://codecov.io/gh/lucasnell/gemino)

# gemino

**An efficient, flexible molecular evolution and sequencing simulator**

## Overview

`gemino` efficiently (i) reads and simulates reference genomes; (ii)
generates variants using summary statistics, phylogenies, Variant Call
Format (VCF) files, and coalescent simulations—the latter of which can
include selection, recombination, and demographic fluctuations; (iii)
simulates sequencing error, mapping qualities, restriction-enzyme
digestion, and variance in coverage among sites; and (iv) writes outputs
to standard file formats. `gemino` can simulate single or paired-ended
reads for WGS on the Illumina platform, and can be extended to simulate
different methods (e.g., original RADseq, double-digest RADseq, and
genotyping-by-sequencing), and sequencing technologies (e.g., Pacific
BioSciences, Oxford Nanopore Technologies).

## Installation

The package is not yet on CRAN, so to install…

``` r
devtools::install_github("lucasnell/gemino")
```

## Usage

``` r
library(gemino)
reference <- create_genome(n_seqs = 10, len_mean = 1000)
phy <- ape::rcoal(5)
mevo_info <- make_mevo(reference, list(model = "JC69", lambda = 0.1))
ref_variants <- create_variants(reference, "phylo", phy, mevo_info)
ref_variants
#>                            << Variants object >>
#> # Variants: 5
#> # Mutations: 21,007
#> 
#>                         << Reference genome info: >>
#> < Set of 10 sequences >
#> # Total size: 10,000 bp
#>   name                          sequence                             length
#> seq0       CCCCAGCGCACCACCGAAATAAAAT...GAATCAAGTATTTATCATCATACTTA      1000
#> seq1       ACAGACCCGTATCTATAATTGAGAG...TTCCAGGAAGCTCTAGAGTTTGACGA      1000
#> seq2       GACAGGTTGAGATCTCACCAGCCGA...ACGGCAGCACCGACAGACTCGGTCAA      1000
#> seq3       GGCGACGCCAGTTGTCTAGACGCAA...CAACGCTGCGCACCATTAGTTCTAAG      1000
#> seq4       ACTTGCTACCCTCTATTGGTTCGTT...GTCAATGGGGGTGAACCGAGTGTCGC      1000
#> seq5       AGCCTGATCTTCGTAATAAGCCACT...TAGTGGGGTTGGGGGGGGTAGAGCCG      1000
#> seq6       AGTGCACGGCTGGTGCTCAACCACG...GGGCGGCAGCTAATGGAATAAGTCCC      1000
#> seq7       TAATTTTGACTGACGAGAGTTGAAC...GTGTACGTTTAGGGGACTTTCTAACG      1000
#> seq8       AATGTGGCTAAAATCAAGCCCGGAA...CAGGTGTACGTGGAAGTCTTTTTATG      1000
#> seq9       CTACTGTTACACCTAACTCGGTATT...TGCGGATCACTCGCGGCAAACAAGAC      1000
```
