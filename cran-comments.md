## Resubmission

This is a resubmission. In this version I have done the following:

* Removed directed quotation marks in the description text

* Added links to Illumina and Pacific Biosciences to the description field of the
  `DESCRIPTION` file

* Added explanation for the abbreviation "PCR" in `DESCRIPTION`. (Note that "FASTA"
  is not an abbreviation: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC280013/>.)

* Added references for sequence simulation methods to the description field in the 
  `DESCRIPTION` file.

* The install section of the `README` file is now based on that from the `ggplot2`
  package so that it conforms to CRAN's policies.

* Remove call to `cat()` in `hts_illumina.R`. All calls to `cat()` are only in print
  methods now.

* Changed `\dontrun{}` to `\donttest{}` in `pacbio.Rd` and `illumina.Rd`.

* In `pacbio.Rd` and `illumina.Rd`, I replaced `\source{}` with `\references{}`.

* Changed title to "A Swift, Versatile Phylogenomic and High-Throughput Sequencing 
  Simulator"



## Test environments
* macOS 10.14.4 (local), R 3.6.0
* ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0
* macOS 10.3.3 (on travis-ci), R-devel, R 3.6.0
* Windows Server 2012 R2 x64 (on appveyor), R 3.6.0
* Windows (win-builder), R-devel, R 3.6.0


## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

```
installed size is 21.0Mb
  sub-directories of 1Mb or more:
    art_profiles   1.5Mb
    libs          18.9Mb
```

The package makes extensive use of compiled code to improve performance.



```
Maintainer: 'Lucas A. Nell <lucas@lucasnell.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Biosciences (11:17, 17:35)
  FASTA (12:47)
  Illumina (10:76, 16:63)
  PacBio (11:30)
  VCF (14:26)
  jackalope (16:6)
  phylogenies (13:65)
```

This is a new R package, but those words are not spelled wrong.


## Downstream dependencies

There are currently no downstream dependencies for this package
