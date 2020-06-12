## Resubmission

This is a resubmission. In this version, I have done the following:

- No errors on installation for R-devel on Windows
- Updated `src/Makevars` to be compatible with using OpenMP on macOS with R >= 4.0.0
- Updated CITATION
- Replaced deprecated `ape::is.binary.tree` with `ape::is.binary.phylo`


## Test environments

* macOS 10.15.4 (local), R 4.0.1
* ubuntu 16.04 (on travis-ci), R-devel, R 4.0.0, R 3.6.3
* macOS 10.13.6 (on travis-ci), R 4.0.1, R 3.6.3
* Windows (win-builder), R-devel, R 4.0.0, R 3.6.3
* Rhub
    - Oracle Solaris 10, R-patched, 32-bit
    - Fedora Linux, R-devel, clang (with valgrind)



## R CMD check results


There were no ERRORs or WARNINGs.


There were 2 NOTEs:

```
installed size is 21.3Mb
  sub-directories of 1Mb or more:
    art_profiles   1.5Mb
    libs          19.1Mb
```

The package makes extensive use of compiled code to improve performance.


```
GNU make is a SystemRequirements.
```

GNU make syntax is required to properly link to the C libraries supplied by the
package `Rhtslib`.




## Downstream dependencies

There are currently no downstream dependencies for this package
