<!-- badges: start -->

[![CRAN
Version](https://img.shields.io/cran/v/bama?style=flat-square&color=blue&label=CRAN)](https://cran.r-project.org/package=bama)

Bayesian Mediation Analysis
===========================

Perform mediation analysis in the presence of high-dimensional mediators
based on the potential outcome framework. Bayesian Mediation Analysis
(BAMA), developed by Song et al (2019) and Song et al (2020), relies on
two Bayesian sparse linear mixed models to simultaneously analyze a
relatively large number of mediators for a continuous exposure and
outcome assuming a small number of mediators are truly active. This
sparsity assumption also allows the extension of univariate mediator
analysis by casting the identification of active mediators as a variable
selection problem and applying Bayesian methods with continuous
shrinkage priors on the effects.

Installation
------------

You can install `bama` via CRAN

    install.packages("bama")

Or devtools

    devtools::install_github("umich-cphds/bama", build_opts = c())

The Github version may contain new features or bug fixes not yet present
on CRAN, so if you are experiencing issues, you may want to try the
Github version of the package.

If you wish to install the package via `devtools`, you will need a C++
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
example (30000, 35000).

    library(bama)

    Y <- bama.data$y
    A <- bama.data$a

    # grab the mediators from the example data.frame
    M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))

    # We just include the intercept term in this example as we have no covariates
    C1 <- matrix(1, 1000, 1)
    C2 <- matrix(1, 1000, 1)
    beta.m  <- rep(0, 100)
    alpha.a <- rep(0, 100)

    out <- bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "BSLMM", seed = 1234,
                burnin = 1000, ndraws = 1100, weights = NULL, inits = NULL, 
                control = list(k = 2, lm0 = 1e-04, lm1 = 1, l = 1))

    # The package includes a function to summarise output from 'bama'
    summary <- summary(out)
    head(summary)

    # Product Threshold Gaussian 
    ptgmod = bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "PTG", seed = 1234,
                  burnin = 1000, ndraws = 1100, weights = NULL, inits = NULL, 
                  control = list(lambda0 = 0.04, lambda1 = 0.2, lambda2 = 0.2))

    mean(ptgmod$beta.a)
    apply(ptgmod$beta.m, 2, mean)
    apply(ptgmod$alpha.a, 2, mean)
    apply(ptgmod$betam_member, 2, mean)
    apply(ptgmod$alphaa_member, 2, mean)

    # Gaussian Mixture Model
    gmmmod = bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "GMM", seed = 1234,
                  burnin = 1000, ndraws = 1100, weights = NULL, inits = NULL, 
                  control = list(phi0 = 0.01, phi1 = 0.01))

    mean(gmmmod$beta.a)
    apply(gmmmod$beta.m, 2, mean)
    apply(gmmmod$alpha.a, 2, mean)

    mean(gmmmod$sigma.sq.a)
    mean(gmmmod$sigma.sq.e)
    mean(gmmmod$sigma.sq.g)

References
==========

Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
dimensional causal mediation effects in omics studies. Biometrics. 2019;
1-11.

Song, Yanyi, Xiang Zhou, Jian Kang, Max T. Aung, Min Zhang, Wei Zhao,
Belinda L. Needham et al. “Bayesian Sparse Mediation Analysis with
Targeted Penalization of Natural Indirect Effects.” arXiv preprint
arXiv:2008.06366 (2020).
