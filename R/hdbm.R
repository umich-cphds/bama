#' @export
hdbm <- function(Y, A, M, C1, C2, beta.m, alpha.a, pi.m, pi.a, burnin = 30000,
                     nsamples = 1000)
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
        stop("beta.m cannot contain missing values")


    if (!is.double(pi.m))
        stop("pi.m must be real number.")
    if ( 0 >= pi.m || pi.m >= 1)
        stop("pi.m must be in (0,1)")
    if (!is.double(pi.a))
        stop("pi.a must be real number.")
    if ( 0 >= pi.a || pi.a >= 1)
        stop("pi.a must be in (0,1)")

  Y <- normalize(Y)
  A <- normalize(A)
  M <- normalize(M)

  run_hdbm_mcmc(Y, A, M, C1, C2, beta.m, alpha.a, pi.m, pi.a, burnin, nsamples)
}

normalize <- function(x)
{
  if (is.vector(x))
    (x - mean(x)) / sd(x)
  else
    apply(x, 2, function(x) (x - mean(x)) / sd(x))
}
