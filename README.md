HDBM

# Install
You can install `hdbm` via devtools:
```
devtools::install_github("umich-cphds/hdbm")
```

# Example

Grab the .RData file from the repo (It will not be downloaded via devtools) and launch R in the directory where the .RData file is. Then,

```
load(".RData")
library(hdbm)

foo = run_hdbm_mcmc(Y = y, A = a, M = m, C1 = c, C2 = c, beta_m = beta.m, alpha_a = alpha.a, pi_m = 0.5, pi_a = 0.5, burnin = 1, nsamples = 1)
```

As it currently stands the mcmc just runs and prints out the values of various paramters for every `burnin` iteration. It technically returns the orginal Y vector, but you can ignore itfor now.

