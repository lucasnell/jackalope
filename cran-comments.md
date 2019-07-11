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

* macOS 10.14.5 (local), R 3.6.0
* ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0
* macOS 10.3.3 (on travis-ci), R 3.6.0
* Windows Server 2012 R2 x64 (on appveyor), R 3.6.0
* Windows (win-builder), R-devel, R 3.6.0


## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

```
installed size is 21.0Mb
  sub-directories of 1Mb or more:
    art_profiles   1.5Mb
    libs          18.9Mb
```

The package makes extensive use of compiled code to improve performance.


## Downstream dependencies

There are currently no downstream dependencies for this package
