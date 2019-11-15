#' @title Estimate parameters of a multivariate elliptical distribution to fit data via Tyler's method
#'
#' @description Estimate parameters of a multivariate elliptical distribution, namely, the mean vector
#' and the covariance matrix, to fit data. Any data sample with NAs will be simply dropped.
#' The algorithm is based on Tyler's method, which normalizes the centered samples to get rid of
#' the shape of the distribution tail. The data is first demeaned (with the geometric mean by default)
#' and normalized. Then the estimation is based on the maximum likelihood estimation (MLE) and the
#' algorithm is obtained from the majorization-minimization (MM) optimization framework.
#' Since Tyler's method can only estimate the covariance matrix up to a scaling factor,
#' a very effective method is employed to recover the scaling factor.
#'
#' @inheritParams fit_mvt
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{cov}: default is the data sample covariance matrix.}}
#'
#' @return A list containing possibly the following elements:
#'         \item{\code{mu}}{Mean vector estimate.}
#'         \item{\code{cov}}{Covariance matrix estimate.}
#'         \item{\code{log_likelihood}}{Value of log-likelihood after converge of the estimation algorithm
#'                                      (only if \code{ftol < Inf}).}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu} and possibly
#'                                       \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (only if \code{return_iterates = TRUE}).}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has reached (\code{FALSE}).}
#'
#' @author Daniel P. Palomar
#'
#' @seealso \code{\link{fit_Cauchy}} and \code{\link{fit_mvt}}
#'
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Tyler’s Scatter Estimator: Existence, Uniqueness, and Algorithms,”
#' IEEE Trans. on Signal Processing, vol. 62, no. 19, pp. 5143-5156, Oct. 2014.
#'
#' @examples
#' library(mvtnorm)       # to generate heavy-tailed data
#' library(fitHeavyTail)
#'
#' X <- rmvt(n = 1000, df = 6)  # generate Student's t data
#' fit_Tyler(X)
#'
#' @importFrom stats cov var
#' @importFrom utils tail
#' @export
fit_Tyler <- function(X, initial = NULL, max_iter = 100, ptol = 1e-3, ftol = Inf, return_iterates = FALSE, verbose = FALSE) {
  ####### error control ########
  X <- try(as.matrix(X), silent = TRUE)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or coercible to a matrix.")
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", 1:ncol(X))
  if (!all(is.na(X) | is.numeric(X))) stop("\"X\" only allows numerical or NA values.")
  if (anyNA(X)) {
    if (verbose) message("X contains NAs, dropping those observations.")
    mask_NA <- apply(X, 1, anyNA)
    X <- X[!mask_NA, , drop = FALSE]
  }
  if (nrow(X) <= 1) stop("Only T=1 sample!!")
  max_iter <- round(max_iter)
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)

  # initialize all parameters
  mu <- if (is.null(initial$mu)) gmean(X, "Gmedian", k = 10) else initial$mu
  if (is.null(initial$cov)) {
    Sigma <- cov(X)
    Sigma <- Sigma/sum(diag(Sigma))
  } else Sigma <- initial$cov
  X_ <- X_ <- X - matrix(mu, T, N, byrow = TRUE)  # demean data
  weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )
  if (ftol < Inf)
    log_likelihood <- (N/2)*sum(log(weights)) - (T/2)*log(det(Sigma))

  # aux function to save iterates
  snapshot <- function() {
    if (ftol < Inf) list(mu = mu, scatter = Sigma, log_likelihood = log_likelihood)
    else list(mu = mu, scatter = Sigma)
  }

  # loop to compute covariance matrix up to a scaling factor with Tyler estimate
  if (return_iterates) iterates_record <- list(snapshot())
  for (iter in 1:max_iter) {
    # record the current status
    Sigma_prev <- Sigma
    if (ftol < Inf) log_likelihood_prev <- log_likelihood

    # Tyler update
    Sigma <- (N/T) * crossprod(sqrt(weights)*X_)  # (N/T) * t(X_) %*% diag(weights) %*% X_
    Sigma <- Sigma/sum(diag(Sigma))
    weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )

    # stopping criterion
    has_param_converged <- all(abs(Sigma - Sigma_prev) <= .5 * ptol * (abs(Sigma) + abs(Sigma_prev)))
    if (ftol < Inf) {
      log_likelihood <- (N/2)*sum(log(weights)) - (T/2)*log(det(Sigma))
      has_fun_converged <- abs(log_likelihood - log_likelihood_prev) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_prev))
    } else has_fun_converged <- TRUE
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    if (has_param_converged && has_fun_converged) break
  }
  if (verbose) message(sprintf("Number of iterations for Tyler estimator = %d\n", iter))

  # finally, recover missing scaling factor
  kappa <- scaling_fitting_ka_with_b(a = diag(Sigma), b = apply(X^2, 2, mean, trim = max(1/T, 0.03)))
  Sigma <- kappa * Sigma

  # return variables
  vars_to_be_returned <- list("mu"  = mu, "cov" = Sigma)
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood <- log_likelihood
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
  }
  vars_to_be_returned$converged <- (iter < max_iter)
  return(vars_to_be_returned)
}




#' @title Estimate parameters of a multivariate elliptical distribution to fit data under a Cauchy distribution
#'
#' @description Estimate parameters of a multivariate elliptical distribution, namely, the mean vector
#' and the covariance matrix, to fit data. Any data sample with NAs will be simply dropped.
#' The estimation is based on the maximum likelihood estimation (MLE) under a Cauchy distribution and
#' the algorithm is obtained from the majorization-minimization (MM) optimization framework.
#' The Cauchy distribution does not have second-order moments and the algorithm actually estimates
#' the scatter matrix. Nevertheless, assuming that the observed data has second-order moments, the
#' covariance matrix is returned by computing the missing scaling factor with a very effective method.
#'
#' @inheritParams fit_mvt
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{cov}: default is the data sample covariance matrix,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix.}}
#'
#' @return A list containing possibly the following elements:
#'         \item{\code{mu}}{Mean vector estimate.}
#'         \item{\code{cov}}{Covariance matrix estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{log_likelihood}}{Value of log-likelihood after converge of the estimation algorithm
#'                                      (only if \code{ftol < Inf}).}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu} and possibly
#'                                       \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (only if \code{return_iterates = TRUE}).}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has reached (\code{FALSE}).}
#'
#' @author Daniel P. Palomar
#'
#' @seealso \code{\link{fit_Tyler}} and \code{\link{fit_mvt}}
#'
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Robust Estimation of Mean and Covariance Matrix Under Heavy-Tailed Distributions,”
#' IEEE Trans. on Signal Processing, vol. 63, no. 12, pp. 3096-3109, June 2015.
#'
#' @examples
#' library(mvtnorm)       # to generate heavy-tailed data
#' library(fitHeavyTail)
#'
#' X <- rmvt(n = 1000, df = 6)  # generate Student's t data
#' fit_Cauchy(X)
#'
#' @importFrom stats cov var
#' @importFrom utils tail
#' @export
fit_Cauchy <- function(X, initial = NULL, max_iter = 100, ptol = 1e-3, ftol = Inf, return_iterates = FALSE, verbose = FALSE) {
  ####### error control ########
  X <- try(as.matrix(X), silent = TRUE)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or coercible to a matrix.")
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", 1:ncol(X))
  if (!all(is.na(X) | is.numeric(X))) stop("\"X\" only allows numerical or NA values.")
  if (anyNA(X)) {
    if (verbose) message("X contains NAs, dropping those observations.")
    mask_NA <- apply(X, 1, anyNA)
    X <- X[!mask_NA, , drop = FALSE]
  }
  if (nrow(X) <= 1) stop("Only T=1 sample!!")
  max_iter <- round(max_iter)
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)

  # initialize all parameters
  mu <- if (is.null(initial$mu)) colMeans(X) else initial$mu
  Sigma <- if (is.null(initial$cov)) cov(X) else initial$cov
  X_ <- X_ <- X - matrix(mu, T, N, byrow = TRUE)  # demean data
  weights <- 1/(1 + rowSums(X_ * (X_ %*% inv(Sigma))))   # 1/( 1 + diag( X_ %*% inv(Sigma) %*% t(X_) ) )
  if (ftol < Inf)
    log_likelihood <- ((N+1)/2)*sum(log(weights)) - (T/2)*log(det(Sigma))

  # aux function to save iterates
  snapshot <- function() {
    if (ftol < Inf) list(mu = mu, scatter = Sigma, log_likelihood = log_likelihood)
    else list(mu = mu, scatter = Sigma)
  }


  # loop
  if (return_iterates) iterates_record <- list(snapshot())
  for (iter in 1:max_iter) {
    # record the current status
    mu_prev <- mu
    Sigma_prev <- Sigma
    if (ftol < Inf) log_likelihood_prev <- log_likelihood

    # update
    mu <- as.vector(weights %*% X)/sum(weights)
    X_ <- X - matrix(mu, T, N, byrow = TRUE)
    beta <- T/(N+1)/sum(weights)  # acceleration
    Sigma <- beta*((N+1)/T) * crossprod(sqrt(weights)*X_)  # (N/T) * t(X_) %*% diag(weights) %*% X_
    weights <- 1/(1 + rowSums(X_ * (X_ %*% inv(Sigma))))   # 1/( 1 + diag( X_ %*% inv(Sigma) %*% t(X_) ) )

    # stopping criterion
    have_params_converged <-
      all(abs(mu - mu_prev)       <= .5 * ptol * (abs(mu) + abs(mu_prev))) &&
      all(abs(Sigma - Sigma_prev) <= .5 * ptol * (abs(Sigma) + abs(Sigma_prev)))
    if (ftol < Inf) {
      log_likelihood <- ((N+1)/2)*sum(log(weights)) - (T/2)*log(det(Sigma))
      has_fun_converged <- abs(log_likelihood - log_likelihood_prev) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_prev))
    } else has_fun_converged <- TRUE
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    if (have_params_converged && has_fun_converged) break
  }
  if (verbose) message(sprintf("Number of iterations for Cauchy estimator = %d\n", iter))

  # finally, recover missing scaling factor for covmat
  kappa <- scaling_fitting_ka_with_b(a = diag(Sigma), b = apply(X^2, 2, mean, trim = max(1/T, 0.03)))

  # return variables
  vars_to_be_returned <- list("mu"  = mu, "cov" = kappa * Sigma, "scatter" = Sigma)
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood <- log_likelihood
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
  }
  vars_to_be_returned$converged <- (iter < max_iter)
  return(vars_to_be_returned)
}





inv <- function(...) solve(...)


#' @importFrom ICSNP spatial.median
Gmedian_of_means <- function(X, k = min(10, ceiling(T/2))) {
  T <- nrow(X)
  N <- ncol(X)
  i1 <- floor(seq(1, T, by = T/k))
  i2 <- c(i1[-1]-1, T)
  mu <- matrix(NA, k, N)
  for (i in 1:k)
    mu[i, ] <- colMeans(X[i1[i]:i2[i], ])
  mu <- ICSNP::spatial.median(X)
  return(mu)
}


#' @importFrom ICSNP spatial.median
gmean <- function(X, method = "mean", k = NULL) {
  switch(method,
         "mean"             = colMeans(X),
         "median"           = apply(X, 2, median),
         "Gmedian"          = ICSNP::spatial.median(X),
         "Gmedian of means" = Gmedian_of_means(X, k),
         stop("Method unknown"))
}


# IRLS method to minimize ||k*a - b||_1 w.r.t. k
scaling_fitting_ka_with_b <- function(a, b, num_iter = 5) {
  for (i in 1:num_iter) {
    w <- if (i == 1) rep(1, length(a))
    else 1 / pmax(1e-5, abs(kappa * a - b))
    kappa <- sum(w*a*b) / sum(w*a^2)
  }
  return(kappa)
}


