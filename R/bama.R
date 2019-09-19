#' High Dimensional Bayesian Mediation
#'
#' `bama` is a Bayesian inference method that uses continuous shrinkage priors
#' for high-dimensional Bayesian mediation analysis, developed by Song et al
#' (2018). \code{bama} provides estimates for the regression coefficients as
#' well as the posterior inclusion probability for ranking mediators.
#'
#' \code{bama} uses two regression models for the two conditional relationships,
#' \eqn{Y | A, M, C} and \eqn{M | A, C}. For the outcome model, \code{bama}
#' uses
#' \deqn{Y = M \beta_M  + A * \beta_A + C* \beta_C + \epsilon_Y}
#' For the mediator model, \code{bama} uses the model
#' \deqn{M = A * \alpha_A + C * \alpha_C + \epsilon_M}
#'
#' For high dimensional tractability, \code{bama} employs continuous Bayesian
#' shrinkage priors to select mediators and makes the two following assumptions:
#' First, it assumes that all the potential mediators contribute small effects
#' in mediating the exposure-outcome relationship. Second, it assumes
#' that only a small proportion of mediators exhibit large effects
#' ("active" mediators). \code{bama} uses a Metropolis-Hastings within Gibbs
#' MCMC to generate posterior samples from the model.
#'
#' @param Y numeric outcome vector
#' @param A numeric exposure vector
#' @param M numeric matrix of mediators of Y and A
#' @param C numeric matrix of extra covariates to include
#' @param beta.m numeric vector of initial beta.m in the outcome model
#' @param alpha.a numeric vector of initial alpha.a in the mediator model
#' @param burnin number of iterations to run the MCMC before sampling
#' @param ndraws number of draws to take from MCMC after the burnin period
#' @return
#' \code{bama} returns a object of type "bama" with 11 elements each of length
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
#' C <- matrix(1, 1000, 1)
#' beta.m  <- rep(0, 100)
#' alpha.a <- rep(0, 100)
#'
#' set.seed(12345)
#' out <- bama(Y, A, M, C, beta.m, alpha.a, burnin = 1000, ndraws = 100)
#'
#' # The package includes a function to summarise output from 'bama'
#' summary <- summary(out)
#' head(summary)
#' @references
#' Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
#' Dimensional Causal Mediation Effects in Omics Studies.
#' bioRxiv \href{https://doi.org/10.1101/467399}{10.1101/467399}
#' @author Alexander Rix
#' @export
bama <- function(Y, A, M, C, beta.m, alpha.a, burnin, ndraws)
{
    if (!is.vector(Y) || !is.numeric(Y))
        stop("'Y' must be a numeric vector.")
    else if (is.integer(Y))
        Y <- as.double(Y)

    if (any(is.na(Y)))
        stop("'Y' must not have missing values.")

    if (!is.vector(A) || !is.numeric(A))
        stop("'A' should be a numeric vector.")
    else if (is.integer(A))
        A <- as.double(A)

    if (any(is.na(A)))
        stop("'A' cannot have missing values.")
    if (length(A) != length(Y))
        stop("Lengths of 'A' and 'Y' do not match.")

    if (!is.matrix(M) || !is.numeric(M))
        stop("'M' must be a numeric matrix.")
    else if (is.integer(M))
        M <- matrix(as.double(M), nrow(M), ncol(M))

    if (any(is.na(M)))
        stop("'M' cannot have missing values.")
    if (nrow(M) != length(Y))
        stop("The number of rows in 'M' does not match the length of 'Y'.")

    if (!is.matrix(C) || !is.numeric(C))
        stop("'C' must be a numeric matrix.")
    else if (is.integer(C))
        C <- matrix(as.double(C), nrow(C), ncol(C))

    if (any(is.na(C)))
        stop("'C' cannot have missing values.")
    if (nrow(C) != length(Y))
        stop("The number of rows in 'M' does not match the length of 'Y'.")

    if (!is.vector(beta.m) || !is.numeric(beta.m))
        stop("'beta.m' must be a numeric vector")
    if (is.integer(beta.m))
        beta.m <- as.double(beta.m)
    if (any(is.na(beta.m)))
        stop("'beta.m' cannot contain missing values.")

    if (!is.vector(alpha.a) || !is.numeric(alpha.a))
        stop("'alpha.a' must be a numeric vector")

    if (is.integer(alpha.a))
        alpha.a <- as.double(alpha.a)
    if (any(is.na(alpha.a)))
        stop("'alpha.a' cannot contain missing values.")

    pi.m <- mean(abs(beta.m)  > 1e-12)
    pi.a <- mean(abs(alpha.a) > 1e-12)

    if (pi.m == 1)
        pi.m <- 0.9
    if (pi.m == 0)
        pi.m <- 0.1
    if (pi.a == 1)
        pi.a <- 0.9
    if (pi.a == 0)
        pi.a <- 0.1

    if (!is.numeric(burnin))
        stop("'burnin' should be a nonnegative integer.")

    if (!is.integer(burnin))
        burnin <- as.integer(burnin)
    if (!is.numeric(ndraws))
        stop("'ndraws' should be a nonnegative integer.")
    if (!is.integer(ndraws))
        ndraws <- as.integer(ndraws)

    bama.out <- run_bama_mcmc(Y, A, M, C, beta.m, alpha.a, pi.m, pi.a,
                                  burnin, ndraws)

    colnames(bama.out$beta.m)  <- colnames(M)
    colnames(bama.out$alpha.a) <- colnames(M)
    colnames(bama.out$r1)      <- colnames(M)
    colnames(bama.out$r3)      <- colnames(M)

    structure(bama.out, class = "bama")
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
summary.bama <- function(object, rank = F, ci = c(0.025, .975), ...)
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
