#
# TO DO:
# 00) Use sum(log(eigen(Sigma)$values)) instead of log(det(Sigma))
# 000) Use obj_value_record <- obj_value_record[1:k]
#      Sigma_diff_record <- Sigma_diff_record[1:k]
# 0) compare with package TTmoment
# 1) include shrinkage to target in most of the functions
# 2) check the 1D case?
# 3) Junyan: extend to NA still with EM
# 4) Junyan: estend to NA with stochastic EM
# 5) Junyan: include shrinkage
#

# Use nu/(nu-2) as variable for checking convergence


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
         "mean" = {
           return(colMeans(X))
         },
         "median" = {
           return(apply(X, 2, median))
         },
         "Gmedian" = {
           return(ICSNP::spatial.median(X))
         },
         "Gmedian of means" = {
           return(Gmedian_of_means(X, k))
         },
         stop("Method unknown")
  )
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
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Tyler’s Scatter Estimator: Existence, Uniqueness, and Algorithms,”
#' IEEE Trans. on Signal Processing, vol. 62, no. 19, pp. 5143-5156, Oct. 2014.
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
momentsTyler <- function(X, verbose = FALSE) {
  max_iter <- 100
  error_th_Sigma <- 1e-3

  #error control
  if (anyNA(X)) stop("This function cannot handle NAs.")
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (T == 1) stop("Only T=1 sample!!")
  if (N == 1) stop("Data is univariate!")

  #first, compute and remove mean
  mu <- gmean(X, "Gmedian", k = 10)
  X_ <- X_ <- X - matrix(mu, T, N, byrow = TRUE)

  #second, compute covariance matrix up to a scaling factor with Tyler estimate
  Sigma <- cov(X_)  #Gaussian initial point
  Sigma <- Sigma/sum(diag(Sigma))

  obj_value_record <-Sigma_diff_record <- rep(NA, max_iter)
  for (k in 1:max_iter) {
    Sigma_prev <- Sigma

    #Tyler update
    weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )
    obj_value_record[k] <- - (N/2)*sum(log(weights)) + (T/2)*log(det(Sigma))
    Sigma <- (N/T) * crossprod( sqrt(weights)*X_ )  # (N/T) * t(X_) %*% diag(weights) %*% X_
    Sigma <- Sigma/sum(diag(Sigma))

    #stopping criterion
    Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
    if (Sigma_diff_record[k] < error_th_Sigma)
      break
  }
  if (verbose)
    cat(sprintf("Number of iterations for Median-Tyler estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(na.omit(Sigma_diff_record), color="green", legend="Sigma") )
  # print( figure(width=700, title="Convergence of objective value", xlab="t", ylab="obj value") %>%
  #          ly_lines(na.omit(obj_value_record), color="blue") )

  #finally, recover missing scaling factor
  sigma2 <- apply(X, 2, var)
  d <- diag(Sigma)
  kappa <- crossprod(sigma2, d)/crossprod(d, d)  # (sigma2 %*% d) / (d %*% d)
  Sigma <- as.numeric(kappa)*Sigma

  return( list(mu = mu, cov = Sigma,
               obj_value_record = na.omit(obj_value_record)) )
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
momentsCauchy <- function(X, verbose = FALSE) {
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
  sigma2 <- apply(X, 2, var)
  d <- diag(Sigma)
  kappa <- crossprod(sigma2, d)/crossprod(d, d)  # (sigma2 %*% d) / (d %*% d)
  Sigma <- as.numeric(kappa)*Sigma

  return( list(mu = mu, cov = Sigma,
               obj_value_record = na.omit(obj_value_record)) )
}



#' Estimation of mean and covariance matrix based on fitting data to the Student-t distribution.
#'
#' The function \code{momentsStudentt} implements a maximum likelihood estimation (MLE) of the mean, covariance matrix,
#' and degrees of freedom nu to fit the data to a Student-t distribution (via the ECME algorithm).
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param nu If passed, then \code{nu} will not be optimized and will be set to that value.
#' @param method Choice of method c("ECM", "ECME") for the maximization of \code{nu}.
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with the following components:
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' \item{\code{nu}}{degrees of freedom}
#' \item{\code{obj_value_record}}{convergence of objective value vs iterations}
#' @author Daniel P. Palomar and Junyan LIU
#' @references
#' Chuanhai Liu and Donald B. Rubin, “ML estimation of the t-distribution using EM and its extensions, ECM and ECME,”
#' Statistica Sinica (5), pp. 19-39, 1995.
#' @examples
#' data(heavy_data)
#' res <- momentsStudentt(heavy_data$X)
#' norm(res$mu - heavy_data$mu, "2")
#' norm(colMeans(heavy_data$X) - heavy_data$mu, "2")
#' norm(res$cov - heavy_data$cov, "F")
#' norm(cov(heavy_data$X) - heavy_data$cov, "F")
#' @export
#' @importFrom stats cov var optimize
#' @importFrom mvtnorm dmvt
momentsStudentt <- function(X, nu = NULL, max_iter = 100, method = "ECM", verbose = FALSE) {
  error_th_nu <- 0.1
  error_th_Sigma <- 1e-3
  error_th_mu <- 1e-3

  #error control
  if (anyNA(X)) stop("This function cannot handle NAs.")
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (T == 1) stop("Only T=1 sample!!")
  if (N == 1) stop("Data is univariate!")

  optimize_nu <- ifelse(is.null(nu), TRUE, FALSE)
  if (!optimize_nu && nu == Inf)
    nu <- 1e15

  #initial point based on sample mean and SCM
  if (optimize_nu)
    nu <- 4
  mu <- colMeans(X)
  Sigma <- (nu-2)/nu*cov(X)  #Sigma is the scatter matrix, not the covariance matrix

  #loop
  nu_record <- obj_value_record <- mu_diff_record <- Sigma_diff_record <- nu_diff_record <- rep(NA, max_iter)
  for (k in 1:max_iter) {
    mu_prev <- mu
    Sigma_prev <- Sigma
    nu_prev <- nu

    # update
    X_ <- X - matrix(mu, T, N, byrow = TRUE)
    tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
    obj_value_record[k] <- sum(mvtnorm::dmvt(X, delta = mu, sigma = Sigma, df = nu, log = TRUE, type = "shifted"))
    # update mu and sigma
    weights <- (nu+N)/(nu+tmp)  #E_tau
    mu <- as.vector(weights %*% X)/sum(weights)
    X_ <- X - matrix(mu, T, N, byrow = TRUE)  #this is slower: sweep(X, 2, FUN = "-", STATS = mu)  #X_ <- X - rep(mu, each = TRUE)  # this is wrong?
    beta <- T/sum(weights)  #acceleration
    X_demeaned_gaussianized <- sqrt(beta) * sqrt(weights) * X_
    Sigma <- (1/T) * crossprod(X_demeaned_gaussianized)  # (1/T) * t(X_) %*% diag(weights) %*% X_
    # update nu
    if (optimize_nu) {
      switch(method,
             "ECM" = {
               # based on minus the Q function of nu
               S <- T*(digamma((N+nu)/2) - log((N+nu)/2)) + sum(log(weights) - weights)  # S is E_log_tau-E_tau
               Q_nu <- function(nu) { - T*(nu/2)*log(nu/2) + T*lgamma(nu/2) - (nu/2)*sum(S) }
               nu <- optimize(Q_nu, interval = c(2 + 1e-16, 100))$minimum
             },
             "ECME" = {
               # based on minus log-likelihood of nu with mu and sigma fixed to mu[k+1] and sigma[k+1]
               tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
               LL_nu <- function(nu) { - sum ( - ((nu+N)/2)*log(nu+tmp) + lgamma( (nu+N)/2 ) - lgamma(nu/2) + (nu/2)*log(nu) ) }
               nu <- optimize(LL_nu, interval = c(2 + 1e-16, 100))$minimum
             },
             stop("Method unknown")
      )
    }

    #stopping criterion
    nu_record[k] <- nu
    nu_diff_record[k] <- abs((nu - nu_prev)/nu_prev)
    #nu_diff_record[k] <- abs((log(nu) - log(nu_prev))/log(nu_prev))
    mu_diff_record[k] <- norm(mu - mu_prev, "2")/norm(mu_prev, "2")
    Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
    if (mu_diff_record[k] < error_th_mu && Sigma_diff_record[k] < error_th_Sigma && nu_diff_record[k] < error_th_nu)
      break
  }
  if (verbose) cat(sprintf("Number of iterations for Student-t estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(na.omit(nu_diff_record), color="black", legend="nu") %>%
  #          ly_lines(na.omit(mu_diff_record), color="blue", legend="mu") %>%
  #          ly_lines(na.omit(Sigma_diff_record), color="green", legend="Sigma") )
  # figure(width=700, title="Convergence of parameters", xlab="t", ylab="nu") %>%
  #   ly_lines(na.omit(nu_record), color="black", legend="nu")

  X_demeaned_gaussianized <- sqrt(nu/(nu - 2)) * X_demeaned_gaussianized
  # note that norm((1/T)*t(X_demeaned_gaussianized) %*% X_demeaned_gaussianized - nu/(nu-2)*Sigma, "F")

  names(mu) <- colnames(X)

  return(list(mu = mu, cov = nu/(nu-2)*Sigma, scatter = Sigma, nu = nu,
              obj_value_record = na.omit(obj_value_record),
              Xcg = X_demeaned_gaussianized))
}
