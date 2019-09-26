# This file is to view the codes spending time
library(covHeavyTail)

set.seed(234)

p <- 100
n <- p * 4
r <- 5
beta <- matrix(rnorm(p*r, 10, 1), p, r)
psi <- rexp(p, 1/10)
sigma <- beta %*% t(beta) + diag(psi, p)
mu <- rnorm(p, 0, 1)
nu <- 7

X <- mvtnorm::rmvt(n = n, sigma = sigma, df = nu, delta = mu)

# a function for adding NAs
addNA <- function(Y, missing_ratio) {
  n <- nrow(Y)
  p <- ncol(Y)
  n_missing <- round(n*missing_ratio)
  missing_mask <- sample(n, n_missing)

  Y_incomp <- Y
  for (i in 1:n_missing) {
    mask_NA <- sample(c(TRUE, FALSE), p, replace = TRUE, prob = c(0.1, 0.9))
    Y_incomp[missing_mask[i], mask_NA] <- NA
  }
  return(Y_incomp)
}
X_wNA <- addNA(X, 0.1)

max_iter <- 100
ptol <- Inf
ftol <- 1e-6

profvis({
  fit_nom <- covTFA(X, ptol = ptol, ftol = ftol, procedure = TRUE)
})

profvis({
  procedure <- TRUE
  initializer <- NULL
  T <- nrow(X)
  p <- ncol(X)
  FA_struct <- ! r == p

  # init all parameters
  alpha <- 1  # an extra variable for PX-EM acceleration
  nu <- ifelse(is.null(initializer$nu), 10, initializer$nu)
  mask_valid <- !is.na(rowSums(X))
  ifelse(is.null(initializer$mu), mu <- colMeans(X[mask_valid, ]), mu <- initializer$mu)
  S <- cov(X[mask_valid, ])
  if (FA_struct) {
    tmp <- eigen(S, symmetric = TRUE)
    ifelse(is.null(initializer$B), B <- tmp$vectors[, 1:r] %*% diag(sqrt(tmp$values[1:r]), r), B <- initializer$B)
    ifelse(is.null(initializer$psi), psi <- diag(S) - diag(B %*% t(B)), psi <- initializer$psi)
    Sigma <- B %*% t(B) + diag(psi, p)
  } else {
    Sigma <- S
  }


  if (ftol < Inf) log_liks <- log_lik <- dmvtWithNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu)


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
      B   <- optB(S = S, r = r, psi_vec = psi)
      psi <- diag(S - B %*% t(B))
      Sigma <- B %*% t(B) + diag(psi, p)
    } else {
      Sigma <- S
    }

    ## ------- check for convergence --------
    if (ptol < Inf) {
      delta_Sigma <- norm(Sigma - Sigma_old, "F") / norm(Sigma_old, "F")
      delta_mu    <- norm(mu - mu_old, "2") / norm(mu_old, "2")
      delta_nu    <- abs(fnu(nu) - fnu(nu_old)) / abs(fnu(nu_old))
      p_convg     <- delta_Sigma < ptol && delta_mu < ptol && delta_nu < ptol
    }

    if (ftol < Inf) {
      log_lik  <- dmvtWithNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu)
      log_liks <- c(log_liks, log_lik)
      delta_loglik <- abs(log_lik_old - log_lik) / abs(log_lik_old)
      f_convg  <- delta_loglik < ftol
    }

    # record the current the variables if required
    if (procedure) proc[[iter + 1]] <- snap()

    if (p_convg && f_convg) break

  }
})
