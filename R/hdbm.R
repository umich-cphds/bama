hdbm <- function(outcome, exposure, mediators, beta.m, alpha.a, pi.m, pi.a,
                   covars.m = NULL, covars.a = NULL, burnin = 30000,
                   n.samples = 1000)
{
  if (!is.vector(outcome))
   stop("outcome should be a numeric vector.")
  if (!is.numeric(outcome))
    stop("outcome must be numeric.")
  else if (is.integer(outcome))
    outcome <- as.double(outcome)
  if (any(is.na(outcome)))
    stop("outcome cannot have missing values.")
  n <- length(outcome)
  outcome <- as.matrix(outcome)

  exposure <- check_input1(exposure, n, "exposure")

  mediators <- check_input2(mediators, n, "mediators")
  q <- ncol(mediators)

  if (!is.null(covars.m)) {
      covars.m  <- check_input2(covars.m, n, "covars.m")
      covars.m  <- normalize(covars.m)
  }
  covars.m  <- cbind(rep(1, n), covars.m)

  if (!is.null(covars.a)) {
      covars.a  <- check_input2(covars.a, n, "covars.a")
      covars.a  <- normalize(covars.a)
  }
  covars.a  <- cbind(rep(1, n), covars.a)

  beta.m  <- check_input1(beta.m, q, "beta.m")
  alpha.a <- check_input1(alpha.a, q, "alpha.a")

  if (!is.double(pi.m))
    stop("pi.m must be real number.")
  if ( 0 > pi.m || pi.m > 1)
    stop("pi.m must be in [0,1]")
  if (!is.double(pi.a))
    stop("pi.a must be real number.")
  if ( 0 > pi.a || pi.a > 1)
    stop("pi.a must be in [0,1]")

  outcome   <- normalize(outcome)
  exposure  <- normalize(exposure)
  mediators <- normalize(mediators)

  # output <- .Call("hdbm_rwrapper", outcome, exposure, mediators, covars.m,
  #                   covars.a, alpha.a, beta.m, pi.m, pi.a, burnin, n.samples)
  names(output) <- c("beta.m", "r1", "alpha.a", "r3", "beta.a", "pi.m", "pi.a",
                       "sigma.m0", "sigma.m1", "sigma.ma0", "sigm.ma1")
  return(output)
}


check_input1 <- function(x, n, name)
{
  if (!is.vector(x))
    stop(sprintf("%s should be a numeric vector.", name))
  if (length(x) != n)
    stop(sprintf("%s does not have the same length as outcome.", name))
  if (!is.numeric(x))
    stop(sprintf("%s must be numeric.", name))
  else if (is.integer(x))
    x <- as.double(x)
  if (any(is.na(x)))
      stop(sprintf("%s cannot have missing values.", name))
  return(as.matrix(x))
}


check_input2 <- function(x, n, name)
{
  if (is.vector(x)) {
    if (!is.numeric(x))
      stop(sprintf("%s must be numeric.", name))
    else if (is.integer(x))
      x <- as.double(x)
    if (any(is.na(x)))
      stop(sprintf("%s cannot have missing values.", name))
    if (length(x) != n)#' @useDynLib hdbm hdbm_rwrapper
      stop(sprintf("%s must have the same length as outcome.", name))
     x <- as.matrix(x)
  }
  else if (is.data.frame(x)) {
    if (nrow(x) != n)
      stop(sprintf("%s must have the same number of rows as the outcome.", name))
    if (sum(sapply(x, is.numeric)) < ncol(x))
        stop(sprintf("%s must be numeric.", name))
    x <- sapply(x, as.double)
    if (any(is.na(x)))
      stop(sprintf("%s cannot contain missing values.", name))
 }
  else if (is.matrix(x)) {
    if (nrow(x) != n)
      stop(sprintf("%s must have the same number of rows as the outcome.", name))
    if (!is.numeric(x))
      stop(sprintf("%s must be numeric.", name))
    x <- matrix(as.double(x), nrow(x), ncol(x))
    if (any(is.na(x)))
      stop(sprintf("%s cannot contain missing values.", name))
  }
  else
    stop(sprintf("%s does not have the correct type.", name))

  return(x)
}

normalize <- function(x)
{
    for (i in 1:ncol(x))
      x[, i] <- (x[, i] - mean(x[, i])) / sd(x[, i])
  x
}
