#
# TO DO:
# 0) compare with package TTmoment
# 1) include shrinkage to target in most of the functions
# 2) check the 1D case?
# 3) Junyan: extend to NA still with EM
# 4) Junyan: estend to NA with stochastic EM
# 5) Junyan: include shrinkage
#


inv <- function(...) solve(...)


Gmedian_of_means <- function(X, k = 10) {
  T <- nrow(X)
  N <- ncol(X)
  i1 <- floor(seq(1, T, by = T/k))
  i2 <- c(i1[-1]-1, T)
  mu <- matrix(NA, k, N)
  for (i in 1:k)
    mu[i, ] <- colMeans(X[i1[i]:i2[i], ])
  mu <- as.vector( Gmedian::Gmedian(mu, init = colMeans(X)) )
  return(mu)
}



#' Estimation of mean (with median) and of covariance matrix via Tyler's estimator.
#'
#' The function \code{momentsTyler} first estimates a robust mean with the median and then (after removing
#' the mean) estimates the covariance matrix based on Tyler's estimator.
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with components
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' @author Daniel P. Palomar
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Tyler’s Scatter Estimator: Existence, Uniqueness, and Algorithms,”
#' IEEE Trans. on Signal Processing, vol. 62, no. 19, pp. 5143-5156, Oct. 2014.
#' @examples
#' #generate heavy-tailed data
#' N <- 40
#' T <- 100
#' mu <- rep(0,N)
#' nv <- 4
#' U <- t(mvtnorm::rmvnorm(n=round(0.7*N), sigma=0.1*diag(N)))
#' R_cov <- U %*% t(U) + diag(N)
#' X <- mvtnorm::rmvt(n=T, delta=mu, sigma=(nv-2)/nv*R_cov, df=nv)
#'
#' #estimate mean and covariance matrix
#' res <- momentsTyler(X)
#' print( res$mu )
#' print( res$cov )
#' @export
#' @importFrom stats cov var
#' @importFrom utils tail
momentsTyler <- function(X, verbose=FALSE) {
  max_iter <- 1000
  error_th_Sigma <- 1e-3

  #error control
  if (anyNA(X)) stop("This function cannot handle NAs.")
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (T == 1) stop("Only T=1 sample!!")
  if (N == 1) stop("Data is univariate!")


  #first, compute and remove mean
  #mu <- colMeans(X)
  #mu <- apply(X, 2, median)
  #mu <- as.vector( Gmedian(X, init=colMeans(X)) )
  mu <- Gmedian_of_means(X, k = 10)
  X_ <- X - rep(mu, each = T)

  #second, compute covariance matrix up to a scaling factor with Tyler estimate
  Sigma <- cov(X_)  #Gaussian initial point
  Sigma <- Sigma/sum(diag(Sigma))
  Sigma_diff_record <- NULL
  obj_value_record <- NULL
  for (k in 1:max_iter) {
    Sigma_prev <- Sigma

    #Tyler update
    weights <- 1/rowSums(X_ * (X_ %*% inv(Sigma)))   # 1/diag( X_ %*% inv(Sigma) %*% t(X_) )
    obj_value_record <- c( obj_value_record, - (N/2)*sum(log(weights)) + (T/2)*log(det(Sigma)) )
    Sigma <- (N/T) * crossprod( sqrt(weights)*X_ )  # (N/T) * t(X_) %*% diag(weights) %*% X_
    Sigma <- Sigma/sum(diag(Sigma))

    #stopping criterion
    Sigma_diff_record <- c(Sigma_diff_record, norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F"))
    if (tail(Sigma_diff_record, 1) < error_th_Sigma)
      break
  }
  if (verbose)
    cat(sprintf("Number of iterations for Median-Tyler estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(Sigma_diff_record, color="green", legend="Sigma") )
  # print( figure(width=700, title="Convergence of objective value", xlab="t", ylab="obj value") %>%
  #          ly_lines(obj_value_record, color="blue") )

  #finally, recover missing scaling factor
  sigma2 <- apply(X, 2, var)
  d <- diag(Sigma)
  kappa <- crossprod(sigma2, d)/crossprod(d, d)  # (sigma2 %*% d) / (d %*% d)
  Sigma <- as.numeric(kappa)*Sigma

  return( list(mu = mu, cov = Sigma) )
}



#' Estimation of mean and covariance matrix based on fitting data to the Cauchy distribution.
#'
#' The function \code{momentsCauchy} implements a maximum likelihood estimation (MLE) of the mean and covariance matrix
#' to fit the data to a Cauchy distribution (Student-t with nv=1).
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with components
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' @author Daniel P. Palomar
#' @references
#' Ying Sun, Prabhu Babu, and Daniel P. Palomar, “Regularized Robust Estimation of Mean and Covariance Matrix Under Heavy-Tailed Distributions,”
#' IEEE Trans. on Signal Processing, vol. 63, no. 12, pp. 3096-3109, June 2015.
#' @examples
#' #generate heavy-tailed data
#' N <- 40
#' T <- 100
#' mu <- rep(0,N)
#' nv <- 4
#' U <- t(mvtnorm::rmvnorm(n=round(0.7*N), sigma=0.1*diag(N)))
#' R_cov <- U %*% t(U) + diag(N)
#' X <- mvtnorm::rmvt(n=T, delta=mu, sigma=(nv-2)/nv*R_cov, df=nv)
#'
#' #estimate mean and covariance matrix
#' res <- momentsCauchy(X)
#' print( res$mu )
#' print( res$cov )
#' @export
#' @importFrom stats cov var
#' @importFrom utils tail
momentsCauchy <- function(X, verbose=FALSE) {
  max_iter <- 1000
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
  Sigma_diff_record <- mu_diff_record <- NULL
  obj_value_record <- NULL
  for( k in 1:max_iter) {
    mu_prev <- mu
    Sigma_prev <- Sigma

    #update
    X_ <- X - rep(mu, each=T)
    weights <- 1/(1+rowSums(X_ * (X_ %*% inv(Sigma))))   # 1/( 1 + diag( X_ %*% inv(Sigma) %*% t(X_) ) )
    obj_value_record <- c( obj_value_record, - ((N+1)/2)*sum(log(weights)) + (T/2)*log(det(Sigma)) )
    mu <- as.vector(weights %*% X)/sum(weights)
    X_ <- X - rep(mu, each=T)  #this could be removed...
    beta <- T/(N+1)/sum(weights)  #acceleration
    Sigma <- beta*((N+1)/T) * crossprod( sqrt(weights)*X_ )  # (N/T) * t(X_) %*% diag(weights) %*% X_

    #stopping criterion
    mu_diff_record <- c(mu_diff_record, norm(mu - mu_prev, "2")/norm(mu_prev, "2"))
    Sigma_diff_record <- c(Sigma_diff_record, norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F"))
    if (tail(mu_diff_record, 1) < error_th_mu &&
        tail(Sigma_diff_record, 1) < error_th_Sigma)
      break
  }
  if (verbose)
    cat(sprintf("Number of iterations for Cauchy estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(mu_diff_record, color="blue", legend="mu") %>%
  #          ly_lines(Sigma_diff_record, color="green", legend="Sigma") )
  # print( figure(width=700, title="Convergence of objective value", xlab="t", ylab="obj value") %>%
  #          ly_lines(obj_value_record, color="blue") )

  #finally, recover missing scaling factor (since we are imposing nu=1 rather than estimating it)
  sigma2 <- apply(X, 2, var)
  d <- diag(Sigma)
  kappa <- crossprod(sigma2, d)/crossprod(d, d)  # (sigma2 %*% d) / (d %*% d)
  Sigma <- as.numeric(kappa)*Sigma

  return( list(mu = mu, cov = Sigma) )
}



#' Estimation of mean and covariance matrix based on fitting data to the Student-t distribution.
#'
#' The function \code{momentsStudentt} implements a maximum likelihood estimation (MLE) of the mean, covariance matrix,
#' and degrees of freedom nv to fit the data to a Student-t distribution (via the ECME algorithm).
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param verbose If \code{TRUE}, info about iterations for convergence will be output.
#' @return A list with components
#' \item{\code{mu}}{mean estimate}
#' \item{\code{cov}}{covariance matrix estimate}
#' \item{\code{nv}}{degrees of freedom}
#' @author Daniel P. Palomar and Junyan LIU
#' @references
#' Chuanhai Liu and Donald B. Rubin, “ML estimation of the t-distribution using EM and its extensions, ECM and ECME,”
#' Statistica Sinica (5), pp. 19-39, 1995.
#' @examples
#' #generate heavy-tailed data
#' N <- 40
#' T <- 100
#' mu <- rep(0,N)
#' nv <- 4
#' U <- t(mvtnorm::rmvnorm(n=round(0.7*N), sigma=0.1*diag(N)))
#' R_cov <- U %*% t(U) + diag(N)
#' X <- mvtnorm::rmvt(n=T, delta=mu, sigma=(nv-2)/nv*R_cov, df=nv)
#'
#' #estimate mean and covariance matrix
#' res <- momentsStudentt(X)
#' print( res$mu )
#' print( res$cov )
#' print( res$nv )
#' @export
#' @importFrom stats cov var optimize
#' @importFrom utils tail
momentsStudentt <- function(X, verbose=FALSE) {
  max_iter <- 1000
  error_th_nv <- 0.1
  error_th_Sigma <- 1e-3
  error_th_mu <- 1e-3

  #error control
  if (anyNA(X)) stop("This function cannot handle NAs.")
  X <- as.matrix(X)
  T <- nrow(X)
  N <- ncol(X)
  if (T == 1) stop("Only T=1 sample!!")
  if (N == 1) stop("Data is univariate!")


  #initial point based on Gaussian
  mu <- colMeans(X)
  Sigma <- cov(X)  #Sigma is the scale matrix, not the covariance matrix
  nv <- 10

  #loop
  mu_diff_record <- Sigma_diff_record <- nv_diff_record <- NULL
  obj_value_record <- NULL
  for (k in 1:max_iter) {
    mu_prev <- mu
    Sigma_prev <- Sigma
    nv_prev <- nv

    #update
    X_ <- X - rep(mu, each = T)
    tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
    obj_value_record <- c( obj_value_record, sum(mvtnorm::dmvt(X, delta = mu, sigma = Sigma, df = nv, log = TRUE, type = "shifted")) )
    #update mu and sigma
    weights <- (nv+N)/(nv+tmp)  #E_tau
    mu <- as.vector(weights %*% X)/sum(weights)
    X_ <- X - rep(mu, each=T)  #this is slower: sweep(X, 2, FUN="-", STATS=mu)
    beta <- T/sum(weights)  #acceleration
    Sigma <- beta*(1/T) * crossprod( sqrt(weights)*X_ )  # (1/T) * t(X_) %*% diag(weights) %*% X_
    #update nu
    # based on minus the Q function of nv
    S <- T*(digamma((N+nv)/2)-log((N+nv)/2)) + sum(log(weights)-weights)  # S is E_log_tau-E_tau
    Q_nv <- function(nv) { - T*(nv/2)*log(nv/2) + T*lgamma(nv/2) - (nv/2)*sum(S) }
    nv <- optimize(Q_nv, interval = c(1e-6, 1e6))$minimum
    # # based on minus log-likelihood of nv with mu and sigma fixed to mu[k+1] and sigma[k+1]
    # tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
    # LL_nv <- function(nv) { - sum ( - ((nv+N)/2)*log(nv+tmp) + lgamma( (nv+N)/2 ) - lgamma(nv/2) + (nv/2)*log(nv) ) }
    # nv <- optimize(LL_nv, interval=c(1e-6, 1e6))$minimum

    #stopping criterion
    nv_diff_record <- c(nv_diff_record, abs((nv - nv_prev)/nv_prev))
    mu_diff_record <- c(mu_diff_record, norm(mu - mu_prev, "2")/norm(mu_prev, "2"))
    Sigma_diff_record <- c(Sigma_diff_record, norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F"))
    if (tail(mu_diff_record, 1) < error_th_mu &&
        tail(Sigma_diff_record, 1) < error_th_Sigma &&
        tail(nv_diff_record, 1) < error_th_nv)
      break
  }
  if (verbose)
    cat(sprintf("Number of iterations for Student-t estimator = %d\n", k))
  # print( figure(width=700, title="Convergence of parameters", xlab="t", ylab="param diff") %>%
  #          ly_lines(nv_diff_record, color="black", legend="nv") %>%
  #          ly_lines(mu_diff_record, color="blue", legend="mu") %>%
  #          ly_lines(Sigma_diff_record, color="green", legend="Sigma") )

  return( list(mu = as.vector(mu), cov = nv/(nv-2)*Sigma, nv = nv) )
}
