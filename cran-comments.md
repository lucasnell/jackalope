
## Test environments

* macOS 10.15 (local), R 3.6.1
* ubuntu 16.04 (on travis-ci), R-devel, R 3.6.1
* macOS 10.3.3 (on travis-ci), R 3.6.1
* Windows Server 2012 R2 x64 (on appveyor), R 3.6.1
* Windows (win-builder), R-devel, R 3.6.1




## R CMD check results


There were no ERRORs or WARNINGs.


There were 2 NOTEs:

```
installed size is 22.0Mb
  sub-directories of 1Mb or more:
    art_profiles   1.5Mb
    libs          19.9Mb
```

The package makes extensive use of compiled code to improve performance.


```
GNU make is a SystemRequirements.
```

GNU make syntax is required to properly link to the C libraries supplied by the
package `Rhtslib`.




## Downstream dependencies

There are currently no downstream dependencies for this package
