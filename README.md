
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.com/lucasnell/jackalope.svg?branch=master)](https://travis-ci.com/lucasnell/jackalope)
[![codecov](https://codecov.io/gh/lucasnell/jackalope/branch/master/graph/badge.svg)](https://codecov.io/gh/lucasnell/jackalope)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/lucasnell/jackalope?branch=master&svg=true)](https://ci.appveyor.com/project/lucasnell/jackalope)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/jackalope)](https://cran.r-project.org/package=jackalope)

# jackalope <img src="man/figures/logo.png" align="right" alt="" width="120" />

## Overview

For studies using high-throughput sequencing (HTS) data, simulations can
be vital for planning sampling design and testing bioinformatic tools.
However, most HTS sequencing tools provide only very simple ways of
adding deviations from a reference genome. For HTS studies that focus on
patterns of genomic variation among individuals, populations, or
species, having a tool that can simulate realistic patterns of molecular
evolution and generate HTS data from those simulations would be quite
useful.

`jackalope` simply and efficiently simulates (i) variants from reference
genomes and (ii) reads from both Illumina and Pacific Biosciences
(PacBio) platforms. It can either read reference genomes from FASTA
files or simulate new ones. Genomic variants can be simulated using
summary statistics, phylogenies, Variant Call Format (VCF) files, and
coalescent simulations—the latter of which can include selection,
recombination, and demographic fluctuations. `jackalope` can simulate
single, paired-end, or mate-pair Illumina reads, as well as reads from
Pacific Biosciences These simulations include sequencing errors, mapping
qualities, multiplexing, and optical/PCR duplicates. All outputs can be
written to standard file formats.

## Installation

``` r
# To install the latest stable version from CRAN:
install.packages("jackalope")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("lucasnell/jackalope")
```

## Usage

Below shows how to simulate a 10kb genome, then create variants from
that genome using a phylogenetic tree:

``` r
library(jackalope)
reference <- create_genome(n_seqs = 10, len_mean = 1000)
tr <- ape::rcoal(5)
ref_variants <- create_variants(reference, vars_phylo(tr), sub_JC69(0.1))
ref_variants
#>                            << Variants object >>
#> # Variants: 5
#> # Mutations: 17,679
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
