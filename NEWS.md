# bama v0.9.1
# Major Changes
* The most notable feature in this release is the renaming of `hdbm` to `bama` to prevent confusion with regards to its relationship with `hdmm`.

## Minor Changes

* Added a summary function to `bama` to calculates the PIP of each mediator, credible intervals, and optionally ranks them.
* Added data quality check to ensure the matrix column norms in M, C1 and C2 are not too small. If they are too small, `NaNs`      could be generated becuase of overflow problems when dividing by the norm.
* Improved documentation

# hdbm v0.9.0

Initial version released
