<!-- README.md is generated from README.Rmd. Please edit that file -->
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
    #>        estimate    ci.lower     ci.upper  pip
    #> m12  0.20147278  0.13564502  0.259956017 1.00
    #> m65 -0.26537533 -0.34471526 -0.197784330 1.00
    #> m89 -0.14310173 -0.22789191 -0.055381855 0.79
    #> m97 -0.02648393 -0.07660828  0.009788986 0.04
    #> m57  0.01766981 -0.02366967  0.055395817 0.03
    #> m86  0.01450784 -0.02348821  0.057174255 0.03

Reference
=========

Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
Dimensional Causal Mediation Effects in Omics Studies. [bioRxiv
467399](https://doi.org/10.1101/467399)
