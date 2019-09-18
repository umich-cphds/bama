<!-- README.md is generated from README.Rmd. Please edit that file -->
bama
====

`bama` is a Bayesian inference method that uses continuous shrinkage
priors for high-dimensional Bayesian mediation analysis, developed by
Song et al (2018). `bama` provides estimates for the regression
coefficients as well as the posterior inclusion probability for ranking
mediators. \#\# Installation You can install `bama` via CRAN

    install.packages("bama")

Or devtools

    devtools::install_github("umich-cphds/bama", build_opts = c())

If you wish to install the package via devtools, you will need a C++
compiler installed. This can be accomplished by installing Rtools on
Windows and Xcode on MacOS.

Example
-------

Taken from the `bama` help file

    library(bama)

    Y <- bama.data$y
    A <- bama.data$a

    # grab the mediators from the example data.frame
    M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))

    # We just include the intercept term in this example.
    C <- matrix(1, nrow(M), 1)
    beta.m  <- rep(0, 100)
    alpha.a <- rep(0, 100)

    set.seed(1245)
    bama.out <- bama(Y, A, M, C, C, beta.m, alpha.a, burnin = 3000, ndraws = 100)

    # Rank mediators and see summary information
    head(summary(bama.out, rank = T))
    #>        estimate     ci.lower    ci.upper  pip
    #> m65 -0.26694986 -0.335969389 -0.18964906 1.00
    #> m12  0.20653140  0.124436993  0.28020939 0.99
    #> m89 -0.15228166 -0.218542178 -0.06832073 0.89
    #> m93  0.04338473 -0.003451208  0.10405256 0.16
    #> m22 -0.03293335 -0.099433981  0.02209072 0.13
    #> m97 -0.03372490 -0.084701429  0.01362444 0.11

Reference
=========

Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
Dimensional Causal Mediation Effects in Omics Studies. [bioRxiv
467399](https://doi.org/10.1101/467399)
