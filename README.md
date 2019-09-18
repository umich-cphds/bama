<!-- README.md is generated from README.Rmd. Please edit that file -->
bama
====

<!-- badges: start -->
<!-- badges: end -->
bama
====

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

    ## Example

    Taken from the `bama` help file

library(bama)

Y &lt;- bama.data*y**A* &lt; −*b**a**m**a*.*d**a**t**a*a

grab the mediators from the example data.frame
==============================================

M &lt;- as.matrix(bama.data\[, paste0("m", 1:100)\], nrow(bama.data))

We just include the intercept term in this example.
===================================================

C &lt;- matrix(1, nrow(M), 1) beta.m &lt;- rep(0, 100) alpha.a &lt;-
rep(0, 100)

set.seed(1245) bama.out &lt;- bama(Y, A, M, C, C, beta.m, alpha.a,
burnin = 3000, ndraws = 100)

Which mediators are active?
===========================

summary(bama.out) \`\`\`

Reference
=========

Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
Dimensional Causal Mediation Effects in Omics Studies. [bioRxiv
467399](https://doi.org/10.1101/467399)
