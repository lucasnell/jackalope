## Resubmission

This is a resubmission. In this version, I have done the following:

- Fixed Rprintf errors inside src/ref_classes.h
- Skip tests that cause Rcpp::stop to be used when using OpenMP



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
Maintainer: ‘Lucas A. Nell <lucnell@gmail.com>’

New maintainer:
  Lucas A. Nell <lucnell@gmail.com>
Old maintainer(s):
  Lucas A. Nell <lucas@lucasnell.com>
```

I have changed my email. I sent a confirmation from the old address
to `CRAN-submissions@R-project.org`.



```
GNU make is a SystemRequirements.
```

GNU make syntax is required to properly link to the C libraries supplied by the
package `Rhtslib`.




## Downstream dependencies

There are currently no downstream dependencies for this package
