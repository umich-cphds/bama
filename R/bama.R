
#' Bayesian Mediation Analysis
#'
#' `bama` is a Bayesian inference method that uses continuous shrinkage priors
#' for high-dimensional Bayesian mediation analysis, developed by Song et al
#' (2019, 2020). \code{bama} provides estimates for the regression coefficients as
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
#' @details NOTE: GMM not currently supported. Instead, use method = 'PTG'.
#'
#' @param Y Length \code{n} numeric outcome vector
#' @param A Length \code{n} numeric exposure vector
#' @param M \code{n x p} numeric matrix of mediators of Y and A
#' @param C1 \code{n x nc1} numeric matrix of extra covariates to include in the
#'     outcome model
#' @param C2 \code{n x nc2} numeric matrix of extra covariates to include in the
#'     mediator model
#' @param method String indicating which method to use. Options are
#' \itemize{
#' \item{"BSLMM"}{ - mixture of two normal components; Song et al. 2019}
#' \item{"PTG"}{ - product threshold Gaussian prior; Song et al. 2020}
#' \item{"GMM"}{ - NOTE: GMM not currently supported. Instead, use method = 'PTG'. four-component Gaussian mixture prior; Song et al. 2020}
#' }
#' @param inits list of initial values for the Gibbs sampler. Options are
#' \itemize{
#' \item{beta.m}{ - Length \code{p} numeric vector of initial \code{beta.m} in the
#'     outcome model. See details for equation}
#' \item{alpha.a}{ - Length \code{p} numeric vector of initial \code{alpha.a} in
#'     the mediator model. See details for equation}
#' }
#' @param burnin number of iterations to run the MCMC before sampling
#' @param seed numeric seed for GIBBS sampler
#' @param ndraws number of draws to take from MCMC (includes burnin draws)
#' @param weights Length \code{n} numeric vector of weights
#' @param control list of Gibbs algorithm control options. These include prior
#' and hyper-prior parameters. Options vary by method selection. If 
#' \code{method = "BSLMM"}
#' \itemize{
#' \item{k}{ - Shape parameter prior for inverse gamma}
#' \item{lm0}{ - Scale parameter prior for inverse gamma for the small normal
#'     components}
#' \item{lm1}{ - Scale parameter prior for inverse gamma for the large normal
#'    components of beta_m}
#' \item{lma1}{ - Scale parameter prior for inverse gamma for the large normal
#'    component of alpha_a}
#' \item{l}{ - Scale parameter prior for the other inverse gamma distributions}
#' }
#' If \code{method = "PTG"}
#' \itemize{
#' \item{lambda0}{ - threshold parameter for product of alpha.a and beta.m effect}
#' \item{lambda1}{ - threshold parameter for beta.m effect}
#' \item{lambda2}{ - threshold parameter for alpha.a effect}
#' \item{ha}{ - inverse gamma shape prior for sigma.sq.a}
#' \item{la}{ - inverse gamma scale prior for sigma.sq.a}
#' \item{h1}{ - inverse gamma shape prior for sigma.sq.e}
#' \item{l1}{ - inverse gamma scale prior for sigma.sq.e}
#' \item{h2}{ - inverse gamma shape prior for sigma.sq.g}
#' \item{l2}{ - inverse gamma scale prior for sigma.sq.g}
#' \item{km}{ - inverse gamma shape prior for tau.sq.b}
#' \item{lm}{ - inverse gamma scale prior for tau.sq.b}
#' \item{kma}{ - inverse gamma shape prior for tau.sq.a}
#' \item{lma}{ - inverse gamma scale prior for tau.sq.a}
#' }
#' If \code{method = "GMM". NOTE: GMM not currently supported. Instead, use method = 'PTG'.}
#' \itemize{
#' \item{phi0}{ - prior beta.m variance}
#' \item{phi1}{ - prior alpha.a variance}
#' \item{a0}{ - prior count of non-zero beta.m and alpha.a effects}
#' \item{a1}{ - prior count of non-zero beta.m and zero alpha.a effects}
#' \item{a2}{ - prior count of zero beta.m and non-zero alpha.a effects}
#' \item{a3}{ - prior count of zero beta.m and zero alpha.a effects}
#' \item{ha}{ - inverse gamma shape prior for sigma.sq.a}
#' \item{la}{ - inverse gamma scale prior for sigma.sq.a}
#' \item{h1}{ - inverse gamma shape prior for sigma.sq.e}
#' \item{l1}{ - inverse gamma scale prior for sigma.sq.e}
#' \item{h2}{ - inverse gamma shape prior for sigma.sq.g}
#' \item{l2}{ - inverse gamma scale prior for sigma.sq.g}
#' }
#' 
#' @return
#' If method = "BSLMM", then \code{bama} returns a object of type "bama" with 12 elements:
#' \describe{
#' \item{beta.m}{\code{ndraws x p} matrix containing outcome model mediator
#'       coefficients.
#' }
#' \item{r1}{\code{ndraws x p} matrix indicating whether or not each beta.m
#'     belongs to the larger normal component (1) or smaller normal
#'     component (0).
#' }
#' \item{alpha.a}{\code{ndraws x p} matrix containing the mediator model
#'     exposure coefficients.
#' }
#' \item{r3}{\code{ndraws x p} matrix indicating whether or not each alpha.a
#'     belongs to the larger normal component (1) or smaller normal component (0).
#' }
#' \item{beta.a}{Vector of length \code{ndraws} containing the beta.a coefficient.}
#' \item{pi.m}{Vector of length \code{ndraws} containing the proportion of
#'     non zero beta.m coefficients.
#' }
#' \item{pi.a}{Vector of length \code{ndraws} containing the proportion of
#'     non zero alpha.a coefficients.
#' }
#'   \item{sigma.m0}{Vector of length \code{ndraws} containing the standard
#'       deviation of the smaller normal component for mediator-outcome
#'       coefficients (beta.m).
#' }
#' \item{sigma.m1}{Vector of length \code{ndraws} containing standard deviation
#'     of the larger normal component for mediator-outcome coefficients (beta.m).
#' }
#' \item{sigma.ma0}{Vector of length \code{ndraws} containing standard
#'     deviation of the smaller normal component for exposure-mediator
#'     coefficients (alpha.a).
#' }
#' \item{sigma.ma1}{Vector of length \code{ndraws} containing standard deviation
#'     of the larger normal component for exposure-mediator coefficients
#'     (alpha.a).
#' }
#' \item{call}{The R call that generated the output.}
#' }
#' 
#' NOTE: GMM not currently supported. Instead, use method = 'PTG'
#' If method = "GMM", then \code{bama} returns a object of type "bama" with:
#' \describe{
#' \item{beta.m}{\code{ndraws x p} matrix containing outcome model mediator
#'       coefficients.}
#' \item{alpha.a}{\code{ndraws x p} matrix containing the mediator model
#'     exposure coefficients.}
#' \item{betam_member}{\code{ndraws x p} matrix of 1's and 0's where
#' item = 1 only if beta.m is non-zero.}
#' \item{alphaa_member}{\code{ndraws x p} matrix of 1's and 0's where
#' item = 1 only if alpha.a is non-zero.}
#' \item{alpha.c}{\code{ndraws x (q2 + p)} matrix containing alpha_c coefficients.
#' Since alpha.c is a matrix of dimension q2 x p, the draws are indexed as
#' alpha_c(w, j) = w * p + j}
#' \item{beta.c}{\code{ndraws x q1} matrix containing beta_c coefficients.
#' Since beta.c is a matrix of dimension q1 x p}
#' \item{beta.a}{Vector of length \code{ndraws} containing the beta.a coefficient.}
#' \item{sigma.sq.a}{Vector of length \code{ndraws} variance of beta.a effect}
#' \item{sigma.sq.e}{Vector of length \code{ndraws} variance of outcome model error}
#' \item{sigma.sq.g}{Vector of length \code{ndraws} variance of mediator model error}
#' }
#' 
#' If method = "PTG", then \code{bama} returns a object of type "bama" with:
#' \describe{
#' \item{beta.m}{\code{ndraws x p} matrix containing outcome model mediator
#'       coefficients.}
#' \item{alpha.a}{\code{ndraws x p} matrix containing the mediator model
#'     exposure coefficients.}
#' \item{alpha.c}{\code{ndraws x (q2 + p)} matrix containing alpha_c coefficients.
#' Since alpha.c is a matrix of dimension q2 x p, the draws are indexed as
#' alpha_c(w, j) = w * p + j}
#' \item{beta.c}{\code{ndraws x q1} matrix containing beta_c coefficients.
#' Since beta.c is a matrix of dimension q1 x p}
#' \item{betam_member}{\code{ndraws x p} matrix of 1's and 0's where
#' item = 1 only if beta.m is non-zero.}
#' \item{alphaa_member}{\code{ndraws x p} matrix of 1's and 0's where
#' item = 1 only if alpha.a is non-zero.}
#' \item{beta.a}{Vector of length \code{ndraws} containing the beta.a coefficient.}
#' \item{sigma.sq.a}{Vector of length \code{ndraws} variance of beta.a effect}
#' \item{sigma.sq.e}{Vector of length \code{ndraws} variance of outcome model error}
#' \item{sigma.sq.g}{Vector of length \code{ndraws} variance of mediator model error}
#' }
#' 
#' 
#' 
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
#' C1 <- matrix(1, 1000, 1)
#' C2 <- matrix(1, 1000, 1)
#' beta.m  <- rep(0, 100)
#' alpha.a <- rep(0, 100)
#'
#' out <- bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "BSLMM", seed = 1234,
#'             burnin = 100, ndraws = 110, weights = NULL, inits = NULL, 
#'             control = list(k = 2, lm0 = 1e-04, lm1 = 1, lma1 = 1, l = 1))
#'
#' # The package includes a function to summarise output from 'bama'
#' summary <- summary(out)
#' head(summary)
#' \donttest{
#' 
#' # Product Threshold Gaussian 
#' ptgmod = bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "PTG", seed = 1234,
#'               burnin = 100, ndraws = 110, weights = NULL, inits = NULL, 
#'               control = list(lambda0 = 0.04, lambda1 = 0.2, lambda2 = 0.2))
#' 
#' mean(ptgmod$beta.a)
#' apply(ptgmod$beta.m, 2, mean)
#' apply(ptgmod$alpha.a, 2, mean)
#' apply(ptgmod$betam_member, 2, mean)
#' apply(ptgmod$alphaa_member, 2, mean)
#' 
# # Gaussian Mixture Model
# gmmmod = bama(Y = Y, A = A, M = M, C1 = C1, C2 = C2, method = "GMM", seed = 1234,
#               burnin = 100, ndraws = 110, weights = NULL, inits = NULL, 
#               control = list(phi0 = 0.01, phi1 = 0.01))
# 
# mean(gmmmod$beta.a)
# apply(gmmmod$beta.m, 2, mean)
# apply(gmmmod$alpha.a, 2, mean)
# 
# mean(gmmmod$sigma.sq.a)
# mean(gmmmod$sigma.sq.e)
# mean(gmmmod$sigma.sq.g)
# 
#' }
#' 
#' @references
#' Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
#' dimensional causal mediation effects in omics studies. Biometrics. 2019;
#' 1-11. \doi{10.1111/biom.13189}
#' 
#' Song, Yanyi, Xiang Zhou, Jian Kang, Max T. Aung, Min Zhang, Wei Zhao, 
#' Belinda L. Needham et al. "Bayesian Sparse Mediation Analysis with 
#' Targeted Penalization of Natural Indirect Effects." 
#' arXiv preprint arXiv:2008.06366 (2020).
#' 
#' @authors Alexander Rix, Michael Kleinsasser
#' @export
bama <- function(Y, A, M, C1, C2, method, burnin, ndraws, weights = NULL, inits = NULL, 
                 control = list(k = 2.0, lm0 = 1e-4, lm1 = 1.0, lma1 = 1.0, l = 1.0, 
                                lambda0 = 0.04, lambda1 = 0.2, lambda2 = 0.2,
                                phi0 = 0.01, phi1 = 0.01, a0 = 0.01 * ncol(M), 
                                a1 = 0.05 * ncol(M), a2 = 0.05 * ncol(M), a3 = 0.89 * ncol(M)), seed = NULL)
{
    call <- match.call()
    
    if (is.null(seed)) 
        seed = floor(runif(1, 0, 1000))

    if (!is.vector(Y) || !is.numeric(Y))
        stop("'Y' must be a numeric vector.")
    if (any(is.na(Y)))
        stop("'Y' must not have missing values.")

    n = length(Y)
    p = ncol(M)
    q1 = ncol(C1)
    q2 = ncol(C2)

    if (!is.vector(A) || !is.numeric(A))
        stop("'A' should be a numeric vector.")
    if (any(is.na(A)))
        stop("'A' cannot have missing values.")
    if (length(A) != n)
        stop("Lengths of 'A' and 'Y' do not match.")

    if (!is.matrix(M) || !is.numeric(M))
        stop("'M' must be a numeric matrix.")
    if (any(is.na(M)))
        stop("'M' cannot have missing values.")
    if (nrow(M) != length(Y))
        stop("The number of rows in 'M' does not match the length of 'Y'.")

    if (!is.matrix(C1) || !is.numeric(C1))
        stop("'C1' must be a numeric matrix.")
    if (any(is.na(C1)))
        stop("'C1' cannot have missing values.")
    if (nrow(C1) != length(Y))
        stop("The number of rows in 'C1' does not match the length of 'Y'.")
    if (!is.matrix(C2) || !is.numeric(C2))
        stop("'C2' must be a numeric matrix.")

    if (any(is.na(C2)))
        stop("'C2' cannot have missing values.")
    if (nrow(C2) != length(Y))
        stop("The number of rows in 'C2' does not match the length of 'Y'.")
    
    if (!("beta.m" %in% names(inits)))
        beta.m = rep(0, ncol(M))
    
    if (!("alpha.a" %in% names(inits)))
        alpha.a = rep(0, ncol(M))

    if (!is.vector(beta.m) || !is.numeric(beta.m))
        stop("'beta.m' must be a numeric vector")
    if (any(is.na(beta.m)))
        stop("'beta.m' cannot contain missing values.")

    if (!is.vector(alpha.a) || !is.numeric(alpha.a))
        stop("'alpha.a' must be a numeric vector")
    if (any(is.na(alpha.a)))
        stop("'alpha.a' cannot contain missing values.")

    if (!is.numeric(burnin) || !is.vector(burnin) || length(burnin) != 1)
        stop("'burnin' should be a nonnegative integer.")

    if (!is.numeric(ndraws) || !is.vector(ndraws) || length(ndraws) != 1)
        stop("'ndraws' should be a nonnegative integer.")

    if (!is.null(weights)) {
        if (!is.numeric(weights) || !is.vector(weights) ||
            length(weights) != n || any(weights < 0))
        {
            stop("'weights' must be a length 'n' nonnegative numeric vector.")
        }

        w  <- sqrt(weights)

        Y  <- w * Y
        A  <- w * A
        M  <- apply(M, 2, function(m) m * w)
        C1 <- apply(C1, 2, function(c1) c1 * w)
        C2 <- apply(C2, 2, function(c2) c2 * w)
    }
    
    bama.out = NULL
    
    if (method == "BSLMM") { 
        
        if (!is.null(seed))
            set.seed(seed)
        
        if (!("k" %in% names(control)))
            k = 2.0
        else
            k = control$k
        if (!("lm0" %in% names(control)))
            lm0 = 1e-4
        else 
            lm0 = control$lm0
        if (!("lm1" %in% names(control)))
            lm1 = 1.0
        else
            lm1 = control$lm1
        if (!("lma1" %in% names(control)))
          lma1 = 1.0
        else
          lma1 = control$lma1 
        if (!("l" %in% names(control)))
            l = 1.0
        else
            l = control$l
        
        if (!is.numeric(k) || !is.vector(k) || length(k) != 1 || k < 0)
            stop("'k' should be a nonnegative number.")
        
        if (!is.numeric(lm0) || !is.vector(lm0) || length(lm0) != 1 || lm0 < 0)
            stop("'lm0' should be a nonnegative number.")
        
        if (!is.numeric(lm1) || !is.vector(lm1) || length(lm1) != 1 || lm1 < 0)
            stop("'lm1' should be a nonnegative number.")
        
        if (!is.numeric(lma1) || !is.vector(lma1) || length(lma1) != 1 || lma1 < 0)
          stop("'lma1' should be a nonnegative number.")
        
        if (!is.numeric(l) || !is.vector(l) || length(l) != 1 || l < 0)
            stop("'l' should be a nonnegative number.")
        
        bama.out <- run_bama_mcmc(Y, A, M, C1, C2, beta.m, alpha.a, burnin, ndraws - burnin,
                                  k, lm0, lm1, lma1, l)
        
        colnames(bama.out$beta.m)  <- colnames(M)
        colnames(bama.out$alpha.a) <- colnames(M)
        colnames(bama.out$r1)      <- colnames(M)
        colnames(bama.out$r3)      <- colnames(M)
        
    } else if (method == "PTG") {
        
        if (!("lambda0" %in% names(control)))
            control$lambda0 = 0.04
        if (!("lambda1" %in% names(control)))
            control$lambda1 = 0.2
        if (!("lambda2" %in% names(control)))
            control$lambda2 = 0.2
        if (!("ha" %in% names(control)))
            control$ha = 2.0
        if (!("la" %in% names(control)))
            control$la = 1.0
        if (!("h1" %in% names(control)))
            control$h1 = 2.0
        if (!("l1" %in% names(control)))
            control$l1 = 1.0
        if (!("h2" %in% names(control)))
            control$h2 = 2.0
        if (!("l2" %in% names(control)))
            control$l2 = 1.0
        if (!("km" %in% names(control)))
            control$km = 1.1
        if (!("lm" %in% names(control)))
            control$lm = 0.09
        if (!("kma" %in% names(control)))
            control$kma = 1.1
        if (!("lma" %in% names(control)))
            control$lma = 0.09
        
        samples_betam = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        samples_alphaa = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        betam_member = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        alphaa_member = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        samples_betaa = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_alphac = matrix(rep(0, (ndraws - burnin) * p * q2), nrow = (ndraws - burnin))
        samples_betac = matrix(rep(0, (ndraws - burnin) * q1), nrow = (ndraws - burnin))
        samples_sigmasqa = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_sigmasqe = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_sigmasqg = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        
        ptg(Y = Y, 
            A = A, 
            M = M, 
            C1 = C1,
            C2 = C2,
            ndraws = ndraws, 
            burnin = burnin, 
            seed = seed,
            lambda0 = control$lambda0,
            lambda1 = control$lambda1,
            lambda2 = control$lambda2,
            ha = control$ha,
            la = control$la,
            h1 = control$h1,
            l1 = control$l1,
            h2 = control$h2,
            l2 = control$l2,
            km = control$km,
            lm = control$lm,
            kma = control$kma,
            lma = control$lma,
            samples_betam = samples_betam, 
            samples_alphaa = samples_alphaa,
            betam_member = betam_member,
            alphaa_member = alphaa_member,
            samples_betaa = samples_betaa,
            samples_alphac = samples_alphac,
            samples_betac = samples_betac,
            samples_sigmasqa = samples_sigmasqa,
            samples_sigmasqe = samples_sigmasqe,
            samples_sigmasqg = samples_sigmasqg)
        
        bama.out$beta.m = samples_betam
        bama.out$alpha.a = samples_alphaa
        bama.out$betam_member = betam_member
        bama.out$alphaa_member = alphaa_member
        bama.out$beta.a = samples_betaa
        bama.out$alpha.c = samples_alphac
        bama.out$beta.c = samples_betac
        bama.out$sigma.sq.a = samples_sigmasqa
        bama.out$sigma.sq.e = samples_sigmasqe
        bama.out$sigma.sq.g = samples_sigmasqg
        
        colnames(bama.out$beta.m)  <- colnames(M)
        colnames(bama.out$alpha.a) <- colnames(M)
        
    } else if (method == "GMM") {
      
      stop("GMM is not currently supported. Instead, use method = 'PTG'.")
        
        if (!("phi0" %in% names(control)))
            control$phi0 = 0.01
        if (!("phi1" %in% names(control)))
            control$phi1 = 0.01
        if (!("a0" %in% names(control)))
            control$a0 = 0.01 * p
        if (!("a1" %in% names(control)))
            control$a1 = 0.05 * p
        if (!("a2" %in% names(control)))
            control$a2 = 0.05 * p
        if (!("a3" %in% names(control)))
            control$a3 = 0.89 * p
        if (!("ha" %in% names(control)))
            control$ha = 2.0
        if (!("la" %in% names(control)))
            control$la = 1.0
        if (!("h1" %in% names(control)))
            control$h1 = 2.0
        if (!("l1" %in% names(control)))
            control$l1 = 1.0
        if (!("h2" %in% names(control)))
            control$h2 = 2.0
        if (!("l2" %in% names(control)))
            control$l2 = 1.0
        
        samples_betam = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        samples_alphaa = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        betam_member = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        alphaa_member = matrix(rep(0, (ndraws - burnin) * p), nrow = (ndraws - burnin))
        samples_betaa = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_alphac = matrix(rep(0, (ndraws - burnin) * p * q2), nrow = (ndraws - burnin))
        samples_betac = matrix(rep(0, (ndraws - burnin) * q1), nrow = (ndraws - burnin))
        samples_sigmasqa = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_sigmasqe = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))
        samples_sigmasqg = matrix(rep(0, (ndraws - burnin) * 1), nrow = (ndraws - burnin))

        # gmm(Y = Y, 
        #     A = A, 
        #     M = M, 
        #     C1 = C1,
        #     C2 = C2,
        #     ndraws = ndraws, 
        #     burnin = burnin, 
        #     seed = seed,
        #     phi0 = control$phi0,
        #     phi1 = control$phi1,
        #     a0 = control$a0,
        #     a1 = control$a1,
        #     a2 = control$a2,
        #     a3 = control$a3,
        #     h1 = control$h1,
        #     l1 = control$l1,
        #     h2 = control$h2,
        #     l2 = control$l2,
        #     ha = control$ha,
        #     la = control$la,
        #     samples_betam = samples_betam, 
        #     samples_alphaa = samples_alphaa,
        #     betam_member = betam_member,
        #     alphaa_member = alphaa_member,
        #     samples_betaa = samples_betaa,
        #     samples_alphac = samples_alphac,
        #     samples_betac = samples_betac,
        #     samples_sigmasqa = samples_sigmasqa,
        #     samples_sigmasqe = samples_sigmasqe,
        #     samples_sigmasqg = samples_sigmasqg)
        
        bama.out$beta.m = samples_betam
        bama.out$alpha.a = samples_alphaa
        bama.out$betam_member = betam_member
        bama.out$alphaa_member = alphaa_member
        bama.out$beta.a = samples_betaa
        bama.out$alpha.c = samples_alphac
        bama.out$beta.c = samples_betac
        bama.out$sigma.sq.a = samples_sigmasqa
        bama.out$sigma.sq.e = samples_sigmasqe
        bama.out$sigma.sq.g = samples_sigmasqg
        
        colnames(bama.out$beta.m)  <- colnames(M)
        colnames(bama.out$alpha.a) <- colnames(M)
        
    } else {
        stop("Selected method is not one of 'BSLMM', 'PTG', or 'GMM'.")
    }
    
    bama.out$call <- call

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
    if (!identical(class(object), "bama"))
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
