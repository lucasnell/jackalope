
# jackalope 1.1.1

* Updated `src/Makevars` to be compatible with using OpenMP on macOS with 
  R >= 4.0.0
* Updated CITATION
* Replaced deprecated `ape::is.binary.tree` with `ape::is.binary.phylo`


# jackalope 1.1.0

* Fixed bug when scaling tree using theta parameter.
* Using "haplotypes" instead of "variants" throughout.
* Updated `R6` class documentation for new `roxygen2` methods.
* Provide method to merge >1 chromosomes, but not the whole genome


# jackalope 1.0.0

* `create_variants` performs better with high mutation rates
* Substitutions now use a random-site, discrete-Gamma model for among-site variability
  (optionally with invariant sites), rather than a fixed-sites model
* No longer depends on `vcfR` package, instead using `htslib` for reading VCF files
* Indel creation now uses "tau-leaping" approximation to the Doob-Gillespie algorithm
* Substitutions are added by calculating transition-probability matrices and sampling
  for substitutions at each site

