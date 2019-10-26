# jackalope 1.0.0

* `create_variants` performs better with high mutation rates
* Substitutions now use a random-site, discrete-Gamma model for among-site variability
  (optionally with invariant sites), rather than a fixed-sites model
* No longer depends on `vcfR` package, instead using `htslib` for reading VCF files
* Indel creation now uses "tau-leaping" approximation to the Doob-Gillespie algorithm
* Substitutions are added by calculating transition-probability matrices and sampling
  for substitutions at each site

