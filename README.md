
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build
Status](https://travis-ci.com/lucasnell/jackal.svg?branch=master)](https://travis-ci.com/lucasnell/jackal)
[![codecov](https://codecov.io/gh/lucasnell/jackal/branch/master/graph/badge.svg)](https://codecov.io/gh/lucasnell/jackal)

# jackal

**An efficient, versatile molecular evolution and sequencing simulator**

## Overview

`jackal` efficiently (i) reads and simulates reference genomes; (ii)
generates variants using summary statistics, phylogenies, Variant Call
Format (VCF) files, and coalescent simulations—the latter of which can
include selection, recombination, and demographic fluctuations; (iii)
simulates sequencing error, mapping qualities, restriction-enzyme
digestion, and variance in coverage among sites; and (iv) writes outputs
to standard file formats. `jackal` can simulate single, paired-end, or
mate-pair Illumina reads, as well as reads from Pacific BioSciences.

## Installation

The package is not yet on CRAN, so to install…

``` r
devtools::install_github("lucasnell/jackal")
```

## Usage

``` r
library(jackal)
reference <- create_genome(n_seqs = 10, len_mean = 1000)
phy <- ape::rcoal(5)
mevo_info <- make_mevo(reference, list(model = "JC69", lambda = 0.1))
ref_variants <- create_variants(reference, "phylo", phy, mevo_info)
ref_variants
#>                            << Variants object >>
#> # Variants: 5
#> # Mutations: 21,301
#> 
#>                         << Reference genome info: >>
#> < Set of 10 sequences >
#> # Total size: 10,000 bp
#>   name                          sequence                             length
#> seq0       CGTACCCCGTCGTGATTCTAGACTC...GTGTGAAAATGGCGATGTTACATGTG      1000
#> seq1       TCTCCTCCTGATACCTGAACTCGCC...ACGATCGCGGGTGTTGCTCCGACCGG      1000
#> seq2       ATCTTAAGGCTCACCACTGAGCCAG...AAATAGCATTAGTCCCGTTCGTAAAA      1000
#> seq3       TCCAACAAGTGGGCCAGCTGCTTTA...ACGGTGTCGCCCAGCGTATTATGATA      1000
#> seq4       ACGATATGCTTTCGCGCCTAGGATC...GAGAGTGTTCTCTTCACATCTGTCCG      1000
#> seq5       GCGTTCTAGTTACGTGTTACCCAGG...GTCATAACATTAGATAATCCATCCGA      1000
#> seq6       CGGCTCTAAACTTACCGGAAGACCA...CGGGACGGCGTGCTTCTAGTAACCGC      1000
#> seq7       CCGCGGGCCAGGACAATCTTAACAT...ACTGAGAGTAATGCGTAGCACTGGCA      1000
#> seq8       ATAACGGAGCTCCACACCGTCCTTG...CGATTGGCTTAAGGGTCCGGAATTAA      1000
#> seq9       ATCCAATTGGACTCGAGGGTGATAG...CACCGGGCACACTAGCGTCCTCGTGA      1000
```
