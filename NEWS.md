# bama v0.9.1
## Major Changes
* The most notable feature in this release is the renaming of `hdbm` to `bama` to prevent confusion with regards to its relationship with `hdmm`.
* bama now only accepts a single extra covariate matrix C, instead of C1 and C2.

## Minor Changes

* Added a summary function to `bama` to calculate the PIP, estimate, and credible interval of each mediator and optionally rank them.
* Added data quality check to ensure the matrix column norms in M, C1 and C2 are not too small. If they are too small, `NaNs`
could be generated because of overflow problems when dividing by the norm.
* Improved documentation
* Added checks to C, the matrix of extra covariates, and alpha.a, the initial
alpha.a values.
# hdbm v0.9.0

Initial version released
