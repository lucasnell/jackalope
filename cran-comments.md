## Resubmission

This is a resubmission. In this version, I have done the following:

- Fixed 3 small bugs to remove some unnecessary copying and compilation warnings
- Updated CITATION


## Test environments

* macOS 11.1 (local), R 4.0.3
* ubuntu 20.04 (GitHub Actions), R-devel, R 4.0.3, R 3.6.3
* macOS 10.15.7 (GitHub Actions), R-devel, R 4.0.3, R 3.6.3
* Windows (win-builder), R-devel, R 4.0.3, R 3.6.3
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
