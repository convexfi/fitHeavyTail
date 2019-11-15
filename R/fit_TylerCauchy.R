#
# TO DO:
# Acceleration?

# 00) Use sum(log(eigen(Sigma)$values)) instead of log(det(Sigma))
# 000) Use obj_value_record <- obj_value_record[1:k]
#      Sigma_diff_record <- Sigma_diff_record[1:k]
# 1) include shrinkage to target in most of the functions
# 2) check the 1D case?



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



#' Estimation of mean (with median) and of covariance matrix via Tyler's estimator.
#'
#' The function \code{momentsTyler} first estimates a robust mean with the median and then (after removing
#' the mean) estimates the covariance matrix based on Tyler's estimator.
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with the following components:
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' \item{\code{obj_value_record}}{convergence of objective value vs iterations}
#' @author Daniel P. Palomar
#'
#' @seealso \code{\link{fit_Cauchy}} and \code{\link{fit_mvt}}
#'
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Tyler’s Scatter Estimator: Existence, Uniqueness, and Algorithms,”
#' IEEE Trans. on Signal Processing, vol. 62, no. 19, pp. 5143-5156, Oct. 2014.
#'
#' @examples
#' data(heavy_data)
#' res <- momentsStudentt(heavy_data$X)
#' norm(res$mu - heavy_data$mu, "2")
#' norm(colMeans(heavy_data$X) - heavy_data$mu, "2")
#' norm(res$cov - heavy_data$cov, "F")
#' norm(cov(heavy_data$X) - heavy_data$cov, "F")
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
  if (ftol < Inf) {
    weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )
    log_likelihood <- - (N/2)*sum(log(weights)) + (T/2)*log(det(Sigma))
  }

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
    weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )
    Sigma <- (N/T) * crossprod( sqrt(weights)*X_ )  # (N/T) * t(X_) %*% diag(weights) %*% X_
    Sigma <- Sigma/sum(diag(Sigma))

    #stopping criterion
    has_param_converged <- all(abs(Sigma - Sigma_prev) <= .5 * ptol * (abs(Sigma_prev) + abs(Sigma)))
    if (ftol < Inf) {
      log_likelihood <- - (N/2)*sum(log(weights)) + (T/2)*log(det(Sigma))
      has_fun_converged <- abs(log_likelihood - log_likelihood_prev) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_prev))
    } else has_fun_converged <- TRUE

    # record the current the variables/loglikelihood if required
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    if (has_param_converged && has_fun_converged) break
  }
  if (verbose) message(sprintf("Number of iterations for Tyler estimator = %d\n", iter))

  # finally, recover missing scaling factor
  kappa <- scaling_fitting_ka_with_b(a = diag(Sigma), b = apply(X^2, 2, mean, trim = max(1/T, 0.03)))
  Sigma <- kappa * Sigma

  # return variables
  vars_to_be_returned <- list("mu"  = mu,
                              "cov" = Sigma)
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood <- log_likelihood
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
  }
  vars_to_be_returned$converged <- (iter < max_iter)
  return(vars_to_be_returned)
}





#' Estimation of mean and covariance matrix based on fitting data to the Cauchy distribution.
#'
#' The function \code{momentsCauchy} implements a maximum likelihood estimation (MLE) of the mean and covariance matrix
#' to fit the data to a Cauchy distribution (Student-t with nu=1).
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with the following components:
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' \item{\code{obj_value_record}}{convergence of objective value vs iterations}
#' @author Daniel P. Palomar
#'
#' @seealso \code{\link{fit_Tyler}} and \code{\link{fit_mvt}}
#'
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Robust Estimation of Mean and Covariance Matrix Under Heavy-Tailed Distributions,”
#' IEEE Trans. on Signal Processing, vol. 63, no. 12, pp. 3096-3109, June 2015.
#' @examples
#' data(heavy_data)
#' res <- momentsStudentt(heavy_data$X)
#' norm(res$mu - heavy_data$mu, "2")
#' norm(colMeans(heavy_data$X) - heavy_data$mu, "2")
#' norm(res$cov - heavy_data$cov, "F")
#' norm(cov(heavy_data$X) - heavy_data$cov, "F")
#' @export
#' @importFrom stats cov var
#' @importFrom utils tail
fit_Cauchy <- function(X, verbose = FALSE) {
  max_iter <- 100
  error_th_mu <- 1e-3
  error_th_Sigma <- 1e-3

  #error control
  if (anyNA(X)) stop("This function cannot handle NAs.")
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (T == 1) stop("Only T=1 sample!!")
  if (N == 1) stop("Data is univariate!")


  #initial point based on Gaussian
  mu <- colMeans(X)
  Sigma <- cov(X)

  #loop
  obj_value_record <- mu_diff_record <- Sigma_diff_record <- rep(NA, max_iter)
  for (k in 1:max_iter) {
    mu_prev <- mu
    Sigma_prev <- Sigma

    #update
    X_ <- X - matrix(mu, T, N, byrow = TRUE)
    weights <- 1/(1+rowSums(X_ * (X_ %*% inv(Sigma))))   # 1/( 1 + diag( X_ %*% inv(Sigma) %*% t(X_) ) )
    obj_value_record[k] <- - ((N+1)/2)*sum(log(weights)) + (T/2)*log(det(Sigma))
    mu <- as.vector(weights %*% X)/sum(weights)
    X_ <- X - matrix(mu, T, N, byrow = TRUE)  #this could be removed...
    beta <- T/(N+1)/sum(weights)  #acceleration
    Sigma <- beta*((N+1)/T) * crossprod( sqrt(weights)*X_ )  # (N/T) * t(X_) %*% diag(weights) %*% X_

    #stopping criterion
    mu_diff_record[k] <- norm(mu - mu_prev, "2")/norm(mu_prev, "2")
    Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
    if (mu_diff_record[k] < error_th_mu &&
        Sigma_diff_record[k] < error_th_Sigma)
      break
  }
  if (verbose)
    cat(sprintf("Number of iterations for Cauchy estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(na.omit(mu_diff_record), color="blue", legend="mu") %>%
  #          ly_lines(na.omit(Sigma_diff_record), color="green", legend="Sigma") )
  # print( figure(width=700, title="Convergence of objective value", xlab="t", ylab="obj value") %>%
  #          ly_lines(na.omit(obj_value_record), color="blue") )

  #finally, recover missing scaling factor (since we are imposing nu=1 rather than estimating it)
  kappa <- scaling_fitting_ka_with_b(a = diag(Sigma), b = apply(X^2, 2, mean, trim = max(1/T, 0.03)))
  Sigma <- kappa * Sigma

  return( list(mu = mu, cov = Sigma,
               obj_value_record = na.omit(obj_value_record)) )
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


