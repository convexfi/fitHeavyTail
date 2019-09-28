# This file is a initial version for functions on estimating the parameters
# when the data is assumed to follow the multivariate Students' t (heavy-tailed) distribution and
# admits a factor model structure

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
#' @references
#' Rui ZHOU, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor analysis parameter estimation,"
#' Lecture Notes in Computer Science, 2020. <doi:10.1109/TSP.2015.2452219>
#'
#' @examples
#' # examples are not yet ready!
#'
#' @export
fit_mvt <- function(X, factors = ncol(X), max_iter = 100, ptol = 1e-3, ftol = Inf, method = "ECM", nu = NULL, initializer = NULL, return_iterates = FALSE) {
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
  optimize_nu <- ifelse(is.null(nu), TRUE, FALSE)
  if (!optimize_nu && nu == Inf) nu <- 1e15  # for numerical stability (for the Gaussian case)
  gamma <- .99
  zeta <- 2e-2

  # initialize all parameters
  alpha <- 1  # an extra variable for PX-EM acceleration
  if (optimize_nu)
    nu <- if (is.null(initializer$nu)) 4
          else initializer$nu
  mu <- if (is.null(initializer$mu)) colMeans(X, na.rm = TRUE)
        else initializer$mu
  SCM <- var(X, na.rm = TRUE)
  if (FA_struct) {  # Sigma is the scale matrix, not the covariance matrix
    SCM_eigen <- eigen(SCM, symmetric = TRUE)
    B <- if (is.null(initializer$B)) SCM_eigen$vectors[, 1:factors] %*% diag(sqrt(SCM_eigen$values[1:factors]), factors)
         else initializer$B
    psi <- if (is.null(initializer$psi)) pmax(0, diag(SCM) - diag(B %*% t(B)))
           else initializer$psi
    Sigma <- (nu-2)/nu * (B %*% t(B) + diag(psi, N))
  } else
    Sigma <- (nu-2)/nu * SCM
  #mask_notNA <- !is.na(rowSums(X))
  if (ftol < Inf) log_likelihood <- ifelse(X_has_NA,
                                           dmvt_withNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu),
                                           sum(mvtnorm::dmvt(X, delta = mu, sigma = Sigma, df = nu, log = TRUE, type = "shifted")))
  snapshot <- function() {
    if (ftol < Inf)
      list(mu = mu, scale = Sigma, nu = nu, log_likelihood = log_likelihood)
    else
      list(mu = mu, scale = Sigma, nu = nu)
  }

  # loop
  if (return_iterates) iterations_record <- list(snapshot())
  for (iter in 1:max_iter) {
    # record the current status
    Sigma_old <- Sigma
    mu_old <- mu
    nu_old <- nu
    if (ftol < Inf) log_likelihood_old <- log_likelihood

    ## -------------- E-step --------------
    if (X_has_NA)
      Q <- Estep(mu, Sigma, psi, nu, X)
    else {
      X_ <- X - matrix(mu, T, N, byrow = TRUE)
      tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
      E_tau <- (nu + N) / (nu + tmp)
      ave_E_tau <- mean(E_tau)
      ave_E_tau_X <- (1/T)*as.vector(E_tau %*% X)
    }

    ## -------------- M-step --------------
    # update mu, alpha, nu
    if (X_has_NA) {
      mu <- Q$ave_E_tau_X / Q$ave_E_tau
      alpha <- Q$ave_E_tau
      Sigma <- Q$ave_E_tau_XX - cbind(mu) %*% rbind(Q$ave_E_tau_X) - cbind(Q$ave_E_tau_X) %*% rbind(mu) + Q$ave_E_tau * cbind(mu) %*% rbind(mu)
      nu  <- optimize_nu(- 1 - Q$ave_E_logtau + Q$ave_E_tau)
    } else {
      mu <- ave_E_tau_X / ave_E_tau
      alpha <- ave_E_tau  # acceleration
      X_ <- X - matrix(mu, T, N, byrow = TRUE)  # this is slower: sweep(X, 2, FUN = "-", STATS = mu)  #X_ <- X - rep(mu, each = TRUE)  # this is wrong?
      ave_E_tau_XX <- (1/T) * crossprod(sqrt(E_tau) * X_)  # (1/T) * t(X_) %*% diag(E_tau) %*% X_
      Sigma <- ave_E_tau_XX / alpha  #TODO{Rui}: this Sigma is divided by alpha, whereas your above on line 102 is not... We need to check

      #TODO{Daniel}: trying to fix the oscillations in the convergence of nu
      # gamma <- .99
      # zeta <- 2e-2
      # gammak <- rep(NA, 100)
      # gammak[1] <- gamma
      # for(k in 2:100)
      #   gammak[k] <- gammak[k-1] * (1 - zeta * gammak[k-1])
      # plot(gammak)
      gamma <- gamma * (1 - zeta * gamma)
      if (optimize_nu)
        nu <- gamma*nu + (1-gamma)*switch(method,
                     "ECM" = {  # based on minus the Q function of nu
                       S <- T*(digamma((N+nu)/2) - log((N+nu)/2)) + sum(log(E_tau) - E_tau)  # S is E_log_tau-E_tau
                       Q_nu <- function(nu) { - T*(nu/2)*log(nu/2) + T*lgamma(nu/2) - (nu/2)*sum(S) }
                       optimize(Q_nu, interval = c(2 + 1e-12, 100))$minimum
                       },
                     "ECME" = {  # based on minus log-likelihood of nu with mu and sigma fixed to mu[k+1] and sigma[k+1]
                       tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
                       LL_nu <- function(nu) { - sum ( - ((nu+N)/2)*log(nu+tmp) + lgamma( (nu+N)/2 ) - lgamma(nu/2) + (nu/2)*log(nu) ) }
                       optimize(LL_nu, interval = c(2 + 1e-12, 100))$minimum
                       },
                     stop("Method unknown."))
    }

    if (X_has_NA) {
      # update B & psi
      S <- Q$ave_E_tau_XX - cbind(mu) %*% rbind(Q$ave_E_tau_X) - cbind(Q$ave_E_tau_X) %*% rbind(mu) + Q$ave_E_tau * cbind(mu) %*% rbind(mu)
      S <- S / alpha
      if (FA_struct) {
        B   <- optB(S = S, factors = factors, psi_vec = psi)
        psi <- pmax(0, diag(S - B %*% t(B)))
        Sigma <- B %*% t(B) + diag(psi, N)
      } else {
        Sigma <- S
      }
    }


    ## -------- stopping criterion --------
    ptol_nu <- 1e-1  #TODO{Daniel}: don't forget to remove this
    have_params_converged <-
      all(abs(mu - mu_old)       <= .5 * ptol * (abs(mu_old) + abs(mu))) &&
      abs(fnu(nu) - fnu(nu_old)) <= .5 * ptol_nu * (abs(fnu(nu_old)) + abs(fnu(nu))) &&
      all(abs(Sigma - Sigma_old) <= .5 * ptol * (abs(Sigma_old) + abs(Sigma)))

    if (ftol < Inf) {
      log_likelihood  <- dmvt_withNA(X = X, delta = mu, sigma = Sigma, df = nu)
      has_fun_converged <- abs(log_likelihood - log_likelihood_old) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_old))
    } else has_fun_converged <- TRUE
    # record the current the variables/loglikelihood if required
    if (return_iterates) iterations_record[[iter + 1]] <- snapshot()
    if (have_params_converged && has_fun_converged) break
  }

  ## -------- return variables --------
  #TODO: colnames, rownames, etc.
  vars_to_be_returned <- list("mu"          = mu,
                              "cov"         = nu/(nu-2) * Sigma,
                              "nu"          = nu,
                              "scale"       = Sigma)
  if (FA_struct) {
    vars_to_be_returned$B   <- B
    vars_to_be_returned$Psi <-  psi
  }
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood <- log_likelihood
  if (return_iterates) {
    names(iterations_record) <- paste("iter", 0:(length(iterations_record)-1))
    vars_to_be_returned$iterations_record <- iterations_record
  }
  return(vars_to_be_returned)
}



fnu <- function(nu) {nu/(nu-2)}


# Expectation step
Estep <- function(mu, full_Sigma, psi, nu, X) {
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
    tmp <- nu + rbind(X_demean[i, mask]) %*% Sigma_inv_obs[[Index]] %*% cbind(X_demean[i, mask])
    E_tau_i <- as.numeric( (nu + sum(mask)) / tmp)
    E_tau <- E_tau + E_tau_i / T

    E_logtau_i <- digamma( (nu+sum(mask))/2 ) - log(tmp/2)
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
    "ave_E_tau"    = E_tau,
    "ave_E_logtau" = as.numeric(E_logtau),
    "ave_E_tau_X"  = E_tau_X,
    "ave_E_tau_XX" = E_tau_XX
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
dmvt_withNA <- function(X, delta, sigma, df) {
  res <- 0
  for (i in 1:nrow(X)) {
    mask <- !is.na(X[i, ])
    tmp <- dmvt(x = X[i, mask], delta = delta[mask], sigma = sigma[mask, mask], df = df)
    res <- res + tmp
  }
  return(res)
}
