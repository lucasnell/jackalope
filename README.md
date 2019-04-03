
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.com/lucasnell/jackalope.svg?branch=master)](https://travis-ci.com/lucasnell/jackalope)
[![codecov](https://codecov.io/gh/lucasnell/jackalope/branch/master/graph/badge.svg)](https://codecov.io/gh/lucasnell/jackalope)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/lucasnell/jackalope?branch=master&svg=true)](https://ci.appveyor.com/project/lucasnell/jackalope)

# jackalope

**An efficient, versatile molecular evolution and sequencing simulator**

## Overview

`jackalope` efficiently (i) reads and simulates reference genomes; (ii)
generates variants using summary statistics, phylogenies, Variant Call
Format (VCF) files, and coalescent simulations—the latter of which can
include selection, recombination, and demographic fluctuations; (iii)
simulates sequencing error, mapping qualities, and optical/PCR
duplicates; and (iv) writes outputs to standard file formats.
`jackalope` can simulate single, paired-end, or mate-pair Illumina
reads, as well as reads from Pacific BioSciences.

## Installation

The package is not yet on CRAN, so to install…

``` r
devtools::install_github("lucasnell/jackalope")
```

## Usage

``` r
library(jackalope)
reference <- create_genome(n_seqs = 10, len_mean = 1000)
phy <- ape::rcoal(5)
mevo_info <- make_mevo(reference, list(model = "JC69", lambda = 0.1))
ref_variants <- create_variants(reference, "phylo", phy, mevo_info)
ref_variants
#>                            << Variants object >>
#> # Variants: 5
#> # Mutations: 17,684
#> 
#>                         << Reference genome info: >>
#> < Set of 10 sequences >
#> # Total size: 10,000 bp
#>   name                          sequence                             length
#> seq0       CTGGCATTGAATCATATGAGGTGGC...GTTGCACGATTGATTAAATTCCTGAA      1000
#> seq1       CACTCCGTCGCACACTAGGTTTCGA...GAGCTCGCGTACATGGAGCATTCTGT      1000
#> seq2       CTTAGCCGGAGCGACTCGGAGCAAC...GCGTAATATGCCAGGTCCCGCGTGGC      1000
#> seq3       CGCCTTCCATTTAGGACTTGTATTG...TAAACTCCATGTGACTGTAATGTCAG      1000
#> seq4       GGGTGATATGGTGTGCATGCTGAAT...AGTCTAGAGTCTCTGGGAGGTCAGGT      1000
#> seq5       TTCGTTGGTGGGTGTCCTATGCTAC...CCCGCCGGTTTGACTTACTCGATTGG      1000
#> seq6       GCATGGACAGATGTGATCTGAGTAT...GACCCCATAAGGCCTGGGACACTGTG      1000
#> seq7       TCGTTTCAACGTCCTTAAGTGTAGT...CTCGTTAGCTCTCCGAGGAGACGAGG      1000
#> seq8       CAGGTAAGTTATCAAAGAACCTTCC...GCATCACCTCGCAAGGAGACTCGTTA      1000
#> seq9       GGTAGTAATTAGGCTTAAAATAGCA...AACAAATGTTCGGCATACGATCTACG      1000
```
