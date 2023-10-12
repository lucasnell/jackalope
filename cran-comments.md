## Resubmission

This is a resubmission. In this version, I have done the following:

- Added required jackalope-package man file
- Updated all docs with roxygen2 v7.2.3
- In CITATION, now using bibentry and for authors, using `c()` on person objects


## Test environments

* macOS 13.4.1 (local), R 4.3.1
* macOS 12.7 (GitHub Actions), R 4.3.1
* ubuntu 22.04.3 (GitHub Actions), R-devel, R 4.3.1, R 4.2.3
* Microsoft Windows Server 2022 10.0.20348 (GitHub Actions), R 4.3.1
* Windows (win-builder), R-devel, R 4.3.1, R 4.2.3
* Rhub
    - Fedora Linux, R-devel, clang (with valgrind)



## R CMD check results


There were no ERRORs or WARNINGs.


There were 2 NOTEs:

```
Specified C++11: please drop specification unless essential
```

The package makes extensive use of C++11 code throughout.


```
GNU make is a SystemRequirements.
```

GNU make syntax is required to properly link to the C libraries supplied by the
package `Rhtslib`.




## Downstream dependencies

There are currently no downstream dependencies for this package
