## Resubmission

This is a resubmission. In this version, I have done the following:

- Fixed a bug related to scaling phylogenetic trees in the `vars_theta` function.
- Updated `R6` class documentation for new `roxygen2` methods.


## Test environments

* macOS 10.15.2 (local), R 3.6.2
* ubuntu 16.04 (on travis-ci), R-devel, R 3.6.2
* macOS 10.3.6 (on travis-ci), R 3.6.2
* Windows (win-builder), R-devel with Rtools40
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
