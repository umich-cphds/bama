#' High Dimensional Bayesian Mediation
#'
#' A Bayesian inference method using continuous shrinkage priors for
#' high-dimensional mediation analysis, developed by Song et al (2018).
#' \code{hdbm} provides estimates for the regression coefficients as well as
#' the posterior inclusion probability for ranking mediators.
#'
#' TODO -- details
#' @param Y numeric response vector.
#' @param A numeric exposure vector.
#' @param M numeric matrix of mediators of Y and A.
#' @param C1 numeric matrix of extra covariates in the response model.
#' @param C2 numeric matrix of extra covariates in the mediator model.
#' @param beta.m numeric vector of initial beta.m in the response model.
#' @param alpha.a numeric vector of initial alpha.a in the mediator model.
#' @param burnin number of iterations to run the MCMC before sampling.
#' @param ndraws number of draws to take from MCMC after the burnin period.
#' @return
#' hdbm returns a list with 11 elements (each of length `nsamples`),
#' sampled from the burned in MCMC:
#' \describe{
#'   \item{beta.m}{Response model mediator coefficients}
#'   \item{r1}{Whether or not each beta.m is active (1) or inactive (0)}
#'   \item{alpha.a}{Rediator model exposure coefficients}
#'   \item{r3}{Whether or not alpha.a is active (1) or inactive (0)}
#'   \item{beta.a}{beta.a coefficient}
#'   \item{pi.m}{Proportion of non zero beta.m coefficients.}
#'   \item{pi.a}{Proportion of non zero alpha.a coefficients.}
#'   \item{sigma.m0}{Standard deviation of the inactive mediators in the
#'       response model.}
#'   \item{sigma.m1}{Standard deviation of the active mediators in the
#'       response model.}
#'   \item{sigma.ma0}{Standard deviation of the inactive mediators in the
#'       mediator model.}
#'   \item{sigma.ma1}{Standard deviation of the active mediators in the
#'       mediator model.}
#' }
#' @examples
#' library(hdbm)
#'
#' Y <- hdbm.data$y
#' A <- hdbm.data$a
#'
#' # grab the mediators from the example data.frame
#' M <- as.matrix(hdbm.data[, paste0("m", 1:100)], nrow(hdbm.data))
#'
#' # We just include the intercept term in this example.
#' C <- matrix(1, nrow(M), 1)
#' beta.m <- rep(0, 100)
#' alpha.a <- rep(0, 100)
#'
#' set.seed(1245)
#' output <- hdbm(Y, A, M, C, C, beta.m, alpha.a, burnin = 3000, nsamples = 100)
#'
#' # Which mediators are active?
#' active <- which(colSums(output$r1 * output$r3) > 50)
#' colnames(M)[active]
#' @references
#' Yanyi Song, Xiang Zhou et al. Bayesian Shrinkage Estimation of High
#' Dimensional Causal Mediation Effects in Omics Studies.
#' [bioRxiv 467399](https://doi.org/10.1101/467399)
#' @export
hdbm <- function(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws)
{
    if (!is.vector(Y) || !is.numeric(Y))
        stop("Y must be a numeric vector.")
    else if (is.integer(Y))
        Y <- as.double(Y)

    if (any(is.na(Y)))
        stop("Y must not have missing values.")

    if (!is.vector(A) || !is.numeric(A))
        stop("A should be a numeric vector.")
    else if (is.integer(A))
        A <- as.double(A)

    if (any(is.na(A)))
        stop("A cannot have missing values.")

    if (length(A) != length(Y))
        stop("Lengths of A and Y do not match.")

    if (!is.matrix(M) || !is.numeric(M))
        stop("M must be a numeric matrix.")
    else if (is.integer(M))
        M <- matrix(as.double(M), nrow(M), ncol(M))

    if (any(is.na(M)))
        stop("M cannot have missing values.")

    if (nrow(M) != length(Y))
        stop("The number of rows in M does not match the length of Y.")

    if (!is.vector(beta.m) || !is.numeric(beta.m))
        stop("beta.m must be a numeric vector")
    if (is.integer(beta.m))
        beta.m <- as.double(beta.m)
    if (any(is.na(beta.m)))
        stop("beta.m cannot contain missing values.")

    pi.m <- mean(abs(beta.m)  > 1e-12)
    pi.a <- mean(abs(alpha.a) > 1e-12)

    if (pi.m == 1 || pi.m == 0)
        pi.m <- 0.5
    if (pi.a == 1 || pi.a == 0)
        pi.a <- 0.5

    if (!is.numeric(burnin))
        stop("burnin should be a nonnegative integer.")
    if (!is.numeric(ndraws))
        stop("ndraws should be a nonnegative integer.")

    # Y <- normalize(Y)
    # A <- normalize(A)
    # M <- normalize(M)

    run_hdbm_mcmc(Y, A, M, C1, C2, beta.m, alpha.a, pi.m, pi.a, burnin, ndraws)
}

# normalize <- function(x)
# {
#   if (is.vector(x))
#     (x - mean(x)) / stats::sd(x)
#   else
#     apply(x, 2, function(x) (x - mean(x)) / sd(x))
# }
