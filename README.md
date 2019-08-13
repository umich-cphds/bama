`hdbm` is currently beta software in the process of being polished, but this
should help you get started.
# Install
You can install `hdbm` via devtools:
```
devtools::install_github("umich-cphds/hdbm")
```

# Example

```
library(hdbm)

output <- hdbm(y, a, m, c, c, beta.m = beta.m, alpha.a = alpha.a,
                 burnin = 10000, nsamples = 100)
```
