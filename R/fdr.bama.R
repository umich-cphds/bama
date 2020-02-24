#' Bayesian Mediation Analysis Controlling For False Discovery
#'
#' \code{fdr.bama} uses the permutation method to estimate the null pip
#' distribution for each mediator and determines a threshold (based off of the
#' \code{dfr} parameter) for significance
#' @param Y Length \code{n} numeric outcome vector
#' @param A Length \code{n} numeric exposure vector
#' @param M \code{n x p} numeric matrix of mediators of Y and A
#' @param beta.m Length \code{p} numeric vector of initial \code{beta.m} in the
#'     outcome model
#' @param alpha.a Length \code{p} numeric vector of initial \code{alpha.a} in
#'     the mediator model
#' @param burnin Number of iterations to run the MCMC before sampling
#' @param ndraws Number of draws to take from MCMC after the burnin period
#' @param C1 \code{n x nc1} numeric matrix of extra covariates to include in the
#'     outcome model
#' @param C2 \code{n x nc2} numeric matrix of extra covariates to include in the
#'     mediator model
#' @param fdr False discovery rate. Default is 0.1
#' @param npermutations The number of permutations to generate while estimating
#'     the null pip distribution. Default is 200
#' @param k Shape parameter prior for inverse gamma. Default is 2.0
#' @param lm0 Scale parameter prior for inverse gamma for the small normal
#'    components. Default is 1e-4
#' @param lm1 Scale parameter prior for inverse gamma for the large normal
#'    components. Default is 1.0
#' @param l Scale parameter prior for the other inverse gamma distributions. Default is 1.0
#' @param mc.cores The number of cores to use while running \code{fdr/.bama}.
#' @param type Type of cluster to make when \code{mc.cores > 1}.
#' @return
#' \code{bama} returns a object of type "fdr.bama" with 11 elements each of length
#' \code{ndraws}), sampled from the burned in MCMC:
#' \describe{
#'   \item{beta.m}{Outcome model mediator coefficients}
#'   \item{r1}{Whether or not each beta.m belongs to the larger normal
#'     component (1) or smaller normal component (0)}
#'   \item{alpha.a}{Mediator model exposure coefficients}
#'   \item{r3}{Whether or not each alpha.a belongs to the larger normal
#'     component (1) or smaller normal component (0)}
#'   \item{beta.a}{beta.a coefficient}
#'   \item{pi.m}{Proportion of non zero beta.m coefficients}
#'   \item{pi.a}{Proportion of non zero alpha.a coefficients}
#'   \item{sigma.m0}{standard deviation of the smaller normal component for
#'     mediator-outcome coefficients (beta.m)}
#'   \item{sigma.m1}{standard deviation of the larger normal component for
#'     mediator-outcome coefficients (beta.m)}
#'   \item{sigma.ma0}{Standard deviation of the smaller normal component for
#'     exposure-mediator coefficients (alpha.a)}
#'   \item{sigma.ma1}{Standard deviation of the larger normal component for
#'     exposure-mediator coefficients (alpha.a)}
#' }
#' @examples
#' library(bama)
#'
#' Y <- bama.data$y
#' A <- bama.data$a
#'
#' # grab the mediators from the example data.frame
#' M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))
#'
#' # We just include the intercept term in this example as we have no covariates
#' # bama defaults C1 and C2 to be a matrix of 1s, so this is purely an
#' # illustration.
#' C <- matrix(1, 1000, 1)
#' beta.m  <- rep(0, 100)
#' alpha.a <- rep(0, 100)
#'
#' set.seed(12345)
#' out <- bama(Y, A, M, C, C, beta.m, alpha.a, burnin = 1000, ndraws = 100)
#'
#' # The package includes a function to summarise output from 'bama'
#' summary <- summary(out)
#' head(summary)
#' @references
#' Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
#' dimensional causal mediation effects in omics studies. Biometrics. 2019;
#' 1-11. [https://doi.org/10.1111/biom.13189]{https://doi.org/10.1111/biom.13189}
#' @author Alexander Rix
#' @export
fdr.bama <- function(Y, A, M, beta.m, alpha.a, burnin, ndraws,
                     C1 = matrix(1, length(Y)), C2 = matrix(1, length(Y)),
                     npermutations = 200, fdr = 0.1, k = 2.0, lm0 = 1e-4,
                     lm1 = 1.0, l = 1.0, mc.cores = 1, type = "PSOCK")
{

    bama.fit <- bama(Y, A, M, beta.m, alpha.a, burnin, ndraws, C1, C2,
                     k, lm0, lm1, l)


    if (npermutations <= 0)
        stop("'npermutations' must be a positive integer.")

    if (fdr <= 0 || fdr >=1)
        stop("'fdr' must be in the interval (0, 1).")

    pi2 <- colMeans(bama.fit$r1 == 0 & bama.fit$r3 == 1)
    pi3 <- colMeans(bama.fit$r1 == 1 & bama.fit$r3 == 0)
    pi4 <- colMeans(bama.fit$r1 == 0 & bama.fit$r3 == 0)

    permute.bama <- function(i) {
        Y.p <- sample(Y)
        bama.r1 <- run_bama_mcmc(Y.p, A, M, C1, C2, beta.m, alpha.a, burnin,
                                 ndraws, k, lm0, lm1, l)

        p2 <- colMeans(bama.r1$r1 == 0 & bama.r1$r3 == 1)


        M.p <- M[sample(nrow(M)),]
        bama.r3 <- run_bama_mcmc(Y, A, M.p, C1, C2, beta.m, alpha.a, burnin,
                                 ndraws, k, lm0, lm1, l)

        p3 <- colMeans(bama.r3$r1 == 1 & bama.r3$r3 == 0)

        Y.p <- sample(Y)
        M.p <- M[sample(nrow(M)),]

        bama.r1.r3 <- run_bama_mcmc(Y.p, A, M.p, C1, C2, beta.m, alpha.a,
                                    burnin, ndraws, k, lm0, lm1, l)

        p4 <- colMeans(bama.r1.r3$r1 == 0 & bama.r1.r3$r3 == 0)

        pi2 / (p2 + p3 + p4) * p2 + pi3 / (p2 + p3 + p4) * p3 +
        pi4 / (p2 + p3 + p4) * p4
    }

    if (mc.cores == 1) {
        fits <- sapply(seq(npermutations), permute.bama)
    }
    else {
        cl <- parallel::makeCluster(mc.cores, type = type)
        parallel::clusterExport(cl, list("Y", "A", "M", "C1", "C2", "beta.m",
                            "alpha.a", "burnin", "ndraws", "k", "lm0", "lm1",
                            "l", "pi2",  "pi3", "pi4"), envir = environment())

        fits <- parallel::parSapply(cl, seq(npermutations), permute.bama)
        parallel::stopCluster(cl)
    }

    fits
}

#' Summarize objects of type "bama"
#'
#' summary.bama summarizes the 'beta.m' estimates from \code{bama} and generates
#' an overall estimate, credible interval, and posterior inclusion probability.
#' @return A data.frame with 4 elements. The beta.m estimates, the estimates'
#'     *credible* interval (which by default is 95\%), and the posterior
#'      inclusion probability (pip) of each 'beta.m'.
#' @param object An object of class "bama".
#' @param rank Whether or not to rank the output by posterior inclusion
#'     probability. Default is TRUE.
#' @param ci The credible interval to calculate. \code{ci} should be a length 2
#'     numeric vector specifying the upper and lower bounds of the CI. By
#'     default, \code{ci = c(0.025, .975)}.
#' @param ... Additional optional arguments to \code{summary}
#' @export
summary.fdr.bama <- function(object, rank = F, ci = c(0.025,, .975), ...)
{
    if (class(object) != "bama")
        stop("'object' is not an bama object.")

    if (!is.logical(rank) || length(rank) != 1)
        stop("'rank' should be a length 1 logical.")

    if (!is.numeric(ci) || length(ci) != 2)
        stop("'ci' should be a length 2 numeric.")
    if (ci[1] >= ci[2])
        stop("'ci[1]' should be less than 'ci[2]'.")

    pip    <- colMeans(object$r1 * object$r3)
    beta.m <- colMeans(object$beta.m)

    credible.l <- apply(object$beta.m, 2, stats::quantile, probs = ci[1])
    credible.h <- apply(object$beta.m, 2, stats::quantile, probs = ci[2])

    out <- data.frame(estimate = beta.m, ci.lower = credible.l,
                          ci.upper = credible.h, pip = pip)
    if (rank)
        out <- out[order(pip, decreasing = T), ]

    out
}

#' Printing bama objects
#'
#' Print a bama object.
#' @param x An object of class 'bama'.
#' @param ... Additional arguments to pass to print.data.frame or summary.bama
#' @export
print.bama <- function(x , ...)
{
    print(summary(x, ...), ...)
}
