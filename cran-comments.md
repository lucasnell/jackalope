## Resubmission

This is a resubmission. In this version I have done the following:

* Explicitly tested for memory-access issues using ASAN, and found none.
* Added `::pcg_detail::` prefix to `extended<...>` in `inst/include/pcg/pcg_random.hpp`
  to avoid Solaris compile error
* Fixed error in `src/io_fasta.cpp` that caused heap buffer overflow as detected in ASAN
* More explicitly created object on stack in `make_ref_genome` function (file
  `seq_classes_access.cpp`) to avoid heap buffer overflow as detected in ASAN
* Initialize C++ pointers to `nullptr`
* Added citation information to `inst/CITATION`


## Test environments

* macOS 10.14.5 (local), R 3.6.1
* ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0
* macOS 10.3.3 (on travis-ci), R 3.6.1
* Windows Server 2012 R2 x64 (on appveyor), R 3.6.1
* Windows (win-builder), R-devel, R 3.6.1


## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

```
Maintainer: 'Lucas A. Nell <lucas@lucasnell.com>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  Biosciences (12:17)
  FASTA (13:47)
  Illumina (10:76, 17:63, 21:16)
  PacBio (12:30, 18:16, 23:5)
  Phylogenomic (3:27)
  SimLoRD (24:5)
  St�cker (24:17)
  VCF (15:26)
  al (22:17, 24:28)
  et (22:14, 24:25)
  phylogenies (14:65)
  polymerase (20:17)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2019-06-29 as memory-access issues
    remained.
```

The package was archived on CRAN due to memory-access issues.
I verified that the current version no longer has these issues by testing on
[a Docker container](https://hub.docker.com/r/rocker/r-devel-ubsan-clang)
that uses the Address Sanitizer and Undefined Behaviour Sanitizer.

Those words are spelled correctly.



```
installed size is 21.0Mb
  sub-directories of 1Mb or more:
    art_profiles   1.5Mb
    libs          18.9Mb
```

The package makes extensive use of compiled code to improve performance.


## Downstream dependencies

There are currently no downstream dependencies for this package
