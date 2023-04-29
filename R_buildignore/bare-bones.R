fit_mvt <- function(X, na_rm = TRUE,
                    nu = c("iterative", "kurtosis", "kurtosis-improved", "MLE-diag", "MLE-diag-resampled", "Hill"),
                    nu_iterative_method = c("POP", "OPP"),
                    initial = NULL,
                    optimize_mu = TRUE,
                    weights = NULL,
                    scale_covmat = TRUE,
                    PX_EM_acceleration = TRUE,
                    factors = ncol(X),
                    max_iter = 100, ptol = 1e-3, ftol = Inf,
                    return_iterates = FALSE, verbose = FALSE) {


  T <- nrow(X)
  N <- ncol(X)

  ############################################
  #           initialize parameters          #
  ############################################
  start_time <- proc.time()[3]
  if (optimize_nu) {  # initial point
    nu <- if (is.null(initial$nu)) nu <- 4  # default initial point
    else initial$nu
  }
  if (!is.numeric(nu)) {
    nu <- switch(nu,
                 "kurtosis"           = nu_from_kurtosis(na.omit(X)),
                 "kurtosis-improved"  = nu_from_kurtosis_improved(na.omit(X)),
                 "MLE-diag"           = nu_mle(na.omit(X), method = "MLE-mv-diagcov"),
                 "MLE-diag-resampled" = nu_mle(na.omit(X), method = "MLE-mv-diagcov-resampled"),
                 "Hill"               = Hill_estimator(na.omit(X)),
                 stop("Method to estimate nu unknown."))
    if (verbose) message(sprintf("Automatically setting nu = %.2f", nu))
  } else if (nu == Inf) nu <- 1e15  # for numerical stability (for the Gaussian case)
  if (optimize_mu) {
    mu <- if (is.null(initial$mu)) colMeans(X, na.rm = TRUE)
    else initial$mu
  } else mu <- rep(0, N)
  if (!is.null(initial$scatter)) Sigma <- initial$scatter
  else {
    Sigma <- if (is.null(initial$cov)) (max(nu, 2.1)-2)/max(nu, 2.1) * var(X, na.rm = TRUE)
    else (max(nu, 2.1)-2)/max(nu, 2.1) * initial$cov
  }




  ###########################
  #           loop          #
  ###########################
  for (iter in 1:max_iter) {
    ## -------------- E-step --------------
      Xc <- X - matrix(mu, T, N, byrow = TRUE)
      r2 <- rowSums(Xc * (Xc %*% solve(Sigma)))  # diag( Xc %*% inv(Sigma) %*% t(Xc) )
      E_tau <- (N + nu) / (nu + r2)  # u_t(r2)
      E_tau <- E_tau * weights
      ave_E_tau <- mean(E_tau)
      ave_E_tau_X <- colMeans(X * E_tau)

    ## -------------- M-step --------------
    # update nu
        nu <- switch(nu_iterative_method,
                     "OPP" = {
                       var_X <- 1/(T-1)*colSums(Xc^2)
                       eta <- sum(var_X)/sum(diag(Sigma))
                       min(getOption("nu_max"), max(getOption("nu_min"), 2*eta/(eta - 1)))
                     },
                     "POP" = {
                       u <- (N + nu)/(nu + r2)
                       r2i <- r2/(1 - r2*u/T)
                       theta <- (1 - N/T) * sum(r2i)/T/N
                       min(getOption("nu_max"), max(getOption("nu_min"), 2*theta/(theta - 1)))
                     })

    # update mu and Sigma (taking into account the new nu)
    E_tau <- (N + nu) / (nu + r2) * weights
    ave_E_tau <- mean(E_tau)
    ave_E_tau_X <- colMeans(X * E_tau)
    if (optimize_mu)
      mu <- ave_E_tau_X / ave_E_tau
    Xc <- X - matrix(mu, T, N, byrow = TRUE)
    ave_E_tau_XX <- (1/T) * crossprod(sqrt(E_tau) * Xc)  # (1/T) * t(Xc) %*% diag(E_tau) %*% Xc
    Sigma <- ave_E_tau_XX

    if (PX_EM_acceleration) {
      alpha <- ave_E_tau
      Sigma <- Sigma / alpha
    }

    ## -------- stopping criterion --------
    ## ...
  }




  # correction factor for scatter and cov matrices
  if (scale_covmat) {
    if (!exists("Xc"))
      Xc <- X - matrix(mu, T, N, byrow = TRUE)
    kappa <- compute_kappa_from_marginals(Xc)
    gamma <- gamma_S_psi1(S = Sigma, T = T, psi1 = 1 + kappa)
    NMSE <- (1/gamma) * (1/T) * (kappa*(2*gamma + N) + gamma + N)
    Sigma <- 1/(NMSE + 1) * Sigma
  }

  # cov matrix
  Sigma_cov <- nu/(nu-2) * Sigma
}
