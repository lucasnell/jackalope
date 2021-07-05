## Resubmission

This is a resubmission. In this version, I have done the following:

- Remove one NULL_ENTRY to support STRICT_R_HEADERS


## Test environments

* macOS 11.4 (local), R 4.1.0
* ubuntu 20.04.2 (GitHub Actions), R-devel, R 4.1.0, R 4.0.5
* macOS 10.15.7 (GitHub Actions), R-devel, R 4.1.0, R 4.0.5
* Windows (win-builder), R-devel, R 4.1.0, R 4.0.5
* Rhub
    - Oracle Solaris 10, R-release, 32-bit
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
