`bama` is a Bayesian inference method that uses continuous shrinkage priors for
Bayesian mediation analysis, developed by Song et al (2018). `bama` provides
estimates for the regression coefficients as well as the posterior inclusion
probability for ranking mediators.

# Install
You can install `bama` via CRAN

# Installation
install.packages("bama")
```
Or devtools
```
devtools::install_github("umich-cphds/bama", build_opts = c())
```

If you wish to install the package via devtools, you will need a C++ compiler
installed. This can be accomplished by installing Rtools on Windows and Xcode
on MacOS.

# Example
Taken from the `bama` help file
```
library(bama)

Y <- bama.data$y
A <- bama.data$a

# grab the mediators from the example data.frame
M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))

# We just include the intercept term in this example.
C <- matrix(1, nrow(M), 1)
beta.m <- rep(0, 100)
alpha.a <- rep(0, 100)

set.seed(1245)
output <- bama(Y, A, M, C, C, beta.m, alpha.a, burnin = 3000, ndraws = 100)

# Which mediators are active?
active <- which(colSums(output$r1 * output$r3) > 50)
colnames(M)[active]
```

# Reference
Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
Dimensional Causal Mediation Effects in Omics Studies.
[bioRxiv 467399](https://doi.org/10.1101/467399)
