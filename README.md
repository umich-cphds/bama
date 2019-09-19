<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![Travis build
status](https://travis-ci.org/umich-cphds/bama.svg?branch=master)](https://travis-ci.org/umich-cphds/bama)
<!-- badges: end -->

Bayesian Mediation Analysis
===========================

`bama` is a Bayesian inference method that uses continuous shrinkage
priors for high-dimensional Bayesian mediation analysis, developed by
Song et al (2018). `bama` provides estimates for the regression
coefficients as well as the posterior inclusion probability for ranking
mediators.

Installation
------------

You can install `bama` via CRAN

    install.packages("bama")

Or devtools

    devtools::install_github("umich-cphds/bama", build_opts = c())

If you wish to install the package via devtools, you will need a C++
compiler installed. This can be accomplished by installing Rtools on
Windows and Xcode on MacOS.

Example
-------

This example is taken from the `bama` help file to help you get started
using the method. Please check the documentation of the function by
typing `?bama::bama`, and the vignette by typing `vingette("bama")` in
R.

`bama` includes an example dataset, `bama.data`. It is a `data.frame`
with a numeric response `y`, numeric exposure `a` and 100 numeric
mediators named `m1, m2, ..., m100`.

We recommend using much larger numbers for `burnin` and `ndraws`, for
example (30000, 1000).

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
    bama.out <- bama(Y, A, M, C, beta.m, alpha.a, burnin = 1000, ndraws = 100)

    # Rank mediators and see summary information
    head(summary(bama.out, rank = T))
    #>         estimate      ci.lower    ci.upper  pip
    #> m12  0.199260031  0.1182283493  0.26774482 1.00
    #> m65 -0.264636309 -0.3447152559 -0.19585711 1.00
    #> m89 -0.142337957 -0.2222596056 -0.05402508 0.77
    #> m86  0.015774996 -0.0235833368  0.05717425 0.04
    #> m93  0.032886535 -0.0009256334  0.07676403 0.04
    #> m48 -0.005354759 -0.0454059793  0.03261661 0.03

Reference
=========

Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
Dimensional Causal Mediation Effects in Omics Studies. [bioRxiv
467399](https://doi.org/10.1101/467399)
