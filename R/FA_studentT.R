# This file is a initial version for functions on estimating the parameters
# when the data is assumed to follow the multivariate Students' t (heavy-tailed) distribution and
# admits a factor model structure

# this method should is available in:
# Rui ZHOU, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor analysis parameter estimation",
# Lecture Notes in Computer Science, 2020

#' @title Robustly estimate parameters of the multivariate Students' t data with (optional) assumption of factor model structure
#'
#' @description Robust paramater estimation of the multivariate Students' t data with (optional) assumption of factor model structure
#'
#' @param X Data matrix
#' @param factors Interger indicating number of factor dimension (default is \code{ncol(X)}, so no factor model assumption).
#' @param max_iter Interger indicating the maximum iterations of estimation method.
#' @param ptol Number (>= 0) indicating the tolerance for parameter changing when judge convergence (default is Inf).
#' @param ftol Number (>= 0) indicating the tolerance for objective changing when judge convergence (default is Inf).
#'             Note: it might be time consuming when use objective changing as a convergence judging criterion, especially when X is high-dimensional.
#' @param initializer A list of initial value of parameters for starting method.
#' @param return_iterates A logical value indicating whether to recode the procedure by iterations.
#'
#' @return The estimated parameters as a list.
#'
#' @author Rui ZHOU and Daniel P. Palomar
#'
#' @examples
#' # examples are not yet ready!
#'
#' @export

covTFA <- function(X, factors = ncol(X), max_iter = 100, ptol = 1e-3, ftol = Inf, initializer = NULL, return_iterates = FALSE) {
  ####### error control ########
  X <- as.matrix(X)
  if (nrow(X) == 1) stop("Only T=1 sample!!")
  if (ncol(X) == 1) stop("Data is univariate!")
  factors <- round(factors)
  max_iter <- round(max_iter)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or can be converted to a matrix.")
  if (factors < 1 || factors > ncol(X)) stop("\"factors\" must satisfy \"1 <= factors <= ncol(X)\"")
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)
  X_has_NA <- anyNA(X)
  FA_struct <- factors != N


  #initial point based on sample mean and SCM
  mu <- colMeans(X)
  Sigma <- (nu-2)/nu * cov(X)  # Sigma is the scale matrix, not the covariance matrix




  # initialize all parameters
  alpha <- 1  # an extra variable for PX-EM acceleration
  nu <- if (is.null(initializer$nu)) 10
        else initializer$nu
  mu <- if (is.null(initializer$mu)) colMeans(X, na.rm = TRUE)
        else initializer$mu
  S <- cov(X, na.rm = TRUE)
  if (FA_struct) {
    S_eigen <- eigen(S, symmetric = TRUE)
    B <- if (is.null(initializer$B)) S_eigen$vectors[, 1:factors] %*% diag(sqrt(S_eigen$values[1:factors]), factors)
         else initializer$B
    psi <- if (is.null(initializer$psi)) pmax(0, diag(S) - diag(B %*% t(B)))
           else initializer$psi
    Sigma <- B %*% t(B) + diag(psi, N)
  } else
    Sigma <- S

  mask_notNA <- !is.na(rowSums(X))
  if (ftol < Inf) log_liks <- log_lik <- dmvt_withNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu)




  # enter loop
  p_convg <- f_convg <- TRUE
  snap <- function() {
    if (ftol < Inf)
      list(mu = mu, Sigma = Sigma, nu = nu, log_lik = log_lik)
    else
      list(mu = mu, Sigma = Sigma, nu = nu)
  }
  proc <- list(snap())

  for (iter in 1:max_iter) {

    # record the current status if necessary
    if (ptol < Inf) {
      Sigma_old <- Sigma
      mu_old <- mu
      nu_old <- nu
    }
    if (ftol < Inf)
      log_lik_old <- log_lik

    ## ------------ E-step -------------------
    expectation <- Estep(mu, Sigma, psi, nu, alpha, X)
    E_tau    <- expectation$E_tau
    E_logtau <- expectation$E_logtau
    E_tau_X  <- expectation$E_tau_X
    E_tau_XX <- expectation$E_tau_XX


    ## ------------ M-step -------------------
    # update mu
    mu <- E_tau_X / E_tau

    # update alpha
    alpha <- E_tau

    # update nu
    tmp <- E_logtau - log(alpha) - (E_tau/alpha)
    nu  <- optimize_nu(-1-tmp)

    # update B & psi
    S <- E_tau_XX - cbind(mu) %*% rbind(E_tau_X) - cbind(E_tau_X) %*% rbind(mu) + E_tau * cbind(mu) %*% rbind(mu)
    if (FA_struct) {
      B   <- optB(S = S, factors = factors, psi_vec = psi)
      psi <- diag(S - B %*% t(B))
      Sigma <- B %*% t(B) + diag(psi, p)
    } else {
      Sigma <- S
    }

    # record the current the variables if required
    if (return_iterates) proc[[iter + 1]] <- snap()


    ## -------- stopping criterion --------
    have_param_converged <-
      all(abs(mu - mu_old)       <= .5 * ptol * (abs(mu_old) + abs(mu))) &&
      abs(fnu(nu) - fnu(nu_old)) <= .5 * ptol * (abs(fnu(nu_old)) + abs(fnu(nu))) &&
      all(abs(Sigma - Sigma_old) <= .5 * ptol * (abs(Sigma_old) + abs(Sigma)))

    if (ftol < Inf) {
      log_lik  <- dmvtWithNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu)
      log_liks <- c(log_liks, log_lik)
      has_fun_converged <- abs(log_lik - log_lik_old) <= .5 * ftol * (abs(log_lik) + abs(log_lik_old))
    } else has_fun_converged <- TRUE

    if (have_param_converged && has_fun_converged) break
  }

  ## -------- return variables --------
  vars_to_be_returned <- list("mu" = mu,
                              "cov" = Sigma / alpha,
                              "nu" = nu,
                              "Sigma_scale" = Sigma)
  if (FA_struct) {
    vars_to_be_returned$B   <- B / sqrt(alpha)
    vars_to_be_returned$Psi <-  psi / alpha
  }
  if (ftol < Inf)
    vars_to_be_returned$log_lik <- log_lik
  if (return_iterates)
    vars_to_be_returned$proc <- proc

  return(vars_to_be_returned)
}



fnu <- function(nu) {nu/(nu-2)}


# Expectation step
Estep <- function(mu, full_Sigma, psi, nu, alpha, X) {
  T <- nrow(X)
  N <- ncol(X)

  # to save memory size and reduce computation time
  # full_Sigma <- B %*% t(B) + diag(psi, N)
  missing_pattern <- unique(!is.na(X))
  Sigma_inv_obs <- list()
  for (i in 1:nrow(missing_pattern)) {
    mask <- missing_pattern[i, ]
    Sigma_inv_obs[[i]] <- solve(full_Sigma[mask, mask])
  }

  # init storage space
  E_tau <- E_logtau <- 0
  E_tau_X <- rep(0, N)
  E_tau_XX <- matrix(0, N, N)

  # calculate...
  X_demean <- X - matrix(mu, T, N, byrow = TRUE)
  for (i in 1:T) {

    # find missing pattern
    mask <- !is.na(X[i, ])
    Index <- indexRowOfMatrix(target_vector = mask, mat = missing_pattern)

    # calculate expectation
    tmp <- nu + alpha * rbind(X_demean[i, mask]) %*% Sigma_inv_obs[[Index]] %*% cbind(X_demean[i, mask])
    E_tau_i <- as.numeric(alpha * (nu + sum(mask)) / tmp)
    E_tau <- E_tau + E_tau_i / T

    E_logtau_i <- log(alpha) + digamma( (nu+sum(mask))/2 ) - log(tmp/2)
    E_logtau <- E_logtau + E_logtau_i / T

    tmp <- X[i, ]
    tmp[!mask] <- mu[!mask] + full_Sigma[!mask, mask] %*% Sigma_inv_obs[[Index]] %*% X_demean[i, mask]
    E_tau_x_i <- E_tau_i * tmp
    E_tau_X <- E_tau_X + E_tau_x_i / T

    E_tau_XX <- E_tau_XX + E_tau_i * cbind(tmp) %*% rbind(tmp)
    E_tau_XX[!mask, !mask] <- E_tau_XX[!mask, !mask] + full_Sigma[!mask, !mask] -
      full_Sigma[!mask, mask] %*% Sigma_inv_obs[[Index]] %*% full_Sigma[mask, !mask]
  }
  E_tau_XX <- E_tau_XX / T

  # correct computation precision error (to be deleted!!!)
  min_eigval <- min(eigen(E_tau_XX)$values)
  diag(E_tau_XX) <- diag(E_tau_XX) - min(0, min_eigval)

  return(list(
    "E_tau"    = E_tau,
    "E_logtau" = E_logtau,
    "E_tau_X"  = E_tau_X,
    "E_tau_XX" = E_tau_XX
  ))
}

# assistant function
indexRowOfMatrix <- function(target_vector, mat) {
  return(which(apply(mat, 1, function(x) all.equal(x, target_vector) == "TRUE")))
}

# solve optimal nu via bisection method:
# log(nu/2) - digamma(nu/2) = y
optimize_nu <- function(y, tol = 1e-4) {
  if (y < 0) stop("y must be positive!")
  L <- 1e-5
  U <- 1000

  while ((log(U)-digamma(U)) > y) U <- U*2
  while ((log(L)-digamma(L)) < y) L <- L*2

  while (1) {
    mid <- (L + U) / 2
    tmp <- log(mid)-digamma(mid) - y

    if (abs(tmp) < tol) break
    if (tmp > 0) L <- mid else U <- mid
  }

  return(2*mid)
}


# a function for computing the optimal B given the psi vector
# see: Lemma 1. in LNCS paper
optB <- function(S, factors, psi_vec) {
  psi_sqrt <- sqrt(psi_vec)
  psi_inv_sqrt <- 1 / psi_sqrt
  tmp <- eigen(S * (psi_inv_sqrt%*%t(psi_inv_sqrt)))
  U <- cbind(tmp$vectors[, 1:factors])
  D <- tmp$values[1:factors]
  Z <- matrix(0, nrow(S), factors)

  for (i in 1:factors) {
    zi <- U[, i]
    Z[, i] <- zi * sqrt(max(1, D[i]) - 1) / norm(zi, "2")
  }

  B <- diag(psi_sqrt) %*% Z;
  return(B)
}

# calculate log-likelihood value with missing data
#' @importFrom mvtnorm dmvt
dmvtWithNA <- function(X, delta, sigma, df) {
  res <- 0
  for (i in 1:nrow(X)) {
    mask <- !is.na(X[i, ])
    tmp <- dmvt(x = X[i, mask], delta = delta[mask], sigma = sigma[mask, mask], df = df)
    res <- res + tmp
  }
  return(res)
}
