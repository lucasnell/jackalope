## Resubmission

This is a resubmission. In this version, I have done the following:

- Removed use of `zlibbioc` package that is no longer available.



## Test environments

* macOS 15.5 (local), R 4.5.1
* macOS 15.6.1 (GitHub Actions), R 4.5.1
* ubuntu 24.04.3 (GitHub Actions), R-devel, R 4.5.1, R 4.4.3
* Microsoft Windows Server 2025 10.0.26100 (GitHub Actions), R 4.5.1
* Windows (win-builder), R-devel, R 4.5.1, R 4.4.3
* Rhub
    - Fedora Linux 38, R-devel, clang (with valgrind)



## R CMD check results


There were no ERRORs or WARNINGs.


There was 2 NOTEs:


```
GNU make is a SystemRequirements.
```

GNU make syntax is required to properly link to the C libraries supplied by the
package `Rhtslib`.


```
installed size is 27.5Mb
sub-directories of 1Mb or more:
  art_profiles   1.5Mb
  libs          25.1Mb
```

The package makes extensive use of compiled code to improve performance.



## Downstream dependencies

There are currently no downstream dependencies for this package
