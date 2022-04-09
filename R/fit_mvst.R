#' @title Estimate parameters of a multivariate (generalized hyperbolic) skewed t distribution to fit data
#'
#' @description Estimate parameters of a multivariate (generalized hyperbolic) skewed Student's t distribution to fit data,
#' namely, the location vector, the scatter matrix, the skewness vector, and the degrees of freedom.
#' The estimation is based on the maximum likelihood estimation (MLE) and the algorithm is
#' obtained from the expectation-maximization (EM) method.
#'
#' @details This function estimates the parameters of a (generalized hyperbolic) multivariate Student's t distribution (\code{mu},
#'          \code{scatter}, \code{gamma} and \code{nu}) to fit the data via the expectation-maximization (EM) algorithm.
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param nu Degrees of freedom of the skewed \eqn{t} distribution (otherwise it will be iteratively estimated).
#' @param gamma Skewness vector of the skewed \eqn{t} distribution (otherwise it will be iteratively estimated).
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{nu}: default is \code{4},}
#'                         \item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{gamma}: default is the sample skewness vector,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         }
#' @param max_iter Integer indicating the maximum number of iterations for the iterative estimation
#'                 method (default is \code{500}).
#' @param ptol Positive number indicating the relative tolerance for the change of the variables
#'             to determine convergence of the iterative method (default is \code{1e-3}).
#' @param ftol Positive number indicating the relative tolerance for the change of the log-likelihood
#'             value to determine convergence of the iterative method (default is \code{Inf}, so it is
#'             not active). Note that using this argument might have a computational cost as a convergence
#'             criterion due to the computation of the log-likelihood (especially when \code{X} is high-dimensional).
#' @param PXEM Logical value indicating whether to use the parameter expansion (PX) EM method to accelerating the convergence.
#' @param return_iterates Logical value indicating whether to record the values of the parameters (and possibly the
#'                        log-likelihood if \code{ftol < Inf}) at each iteration (default is \code{FALSE}).
#' @param verbose Logical value indicating whether to allow the function to print messages (default is \code{FALSE}).
#'
#' @return A list containing (possibly) the following elements:
#'         \item{\code{mu}}{Location vector estimate (not the mean).}
#'         \item{\code{gamma}}{Skewness vector estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{mean}}{Mean vector estimate: \preformatted{  mean = mu + nu/(nu-2) * gamma}}
#'         \item{\code{cov}}{Covariance matrix estimate: \preformatted{  cov = nu/(nu-2) * scatter + 2*nu^2 / (nu-2)^2 / (nu-4) * gamma*gamma'}}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has been reached (\code{FALSE}).}
#'         \item{\code{num_iterations}}{Number of iterations executed.}
#'         \item{\code{cpu_time}}{Elapsed overall CPU time.}
#'         \item{\code{log_likelihood_vs_iterations}}{Value of log-likelihood over the iterations (if \code{ftol < Inf}).}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu}, \code{scatter}, \code{nu},
#'                                       and possibly \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (if \code{return_iterates = TRUE}).}
#'         \item{\code{cpu_time_at_iter}}{Elapsed CPU time at each iteration (if \code{return_iterates = TRUE}).}
#'
#'
#' @author Rui Zhou, Xiwen Wang, and Daniel P. Palomar
#'
#' @seealso \code{\link{fit_mvt}}
#'
#' @references
#' Aas Kjersti and Ingrid Hobæk Haff. "The generalized hyperbolic skew Student’s t-distribution,"
#' Journal of financial econometrics, pp. 275-309, 2006.
#'
#' @examples
#' library(mvtnorm)       # to generate heavy-tailed data
#' library(fitHeavyTail)
#'
#' # parameter setting
#' N <- 5
#' T <- 200
#' nu <- 6
#' mu <- rnorm(N)
#' scatter <- diag(N)
#' gamma <- rnorm(N)   # skewness vector
#'
#' # generate GH Skew t data
#' taus <- rgamma(n = T, shape = nu/2, rate = nu/2)
#' X <- matrix(data = mu, nrow = T, ncol = N, byrow = TRUE) +
#'      matrix(data = gamma, nrow = T, ncol = N, byrow = TRUE) / taus +
#'      rmvnorm(n = T, mean = rep(0, N), sigma = scatter) / sqrt(taus)
#'
#' # fit skew t model
#' fit_mvst(X)
#'
#' # setting lower limit for nu (e.g., to guarantee existence of co-skewness and co-kurtosis matrices)
#' options(nu_min = 8.01)
#' fit_mvst(X)
#'
#' @importFrom stats optimize
#' @export
fit_mvst <- function(X,
                     nu = NULL, gamma = NULL, initial = NULL,
                     max_iter = 500, ptol = 1e-3, ftol = Inf,
                     PXEM = TRUE, return_iterates = FALSE, verbose = FALSE) {
  ####### error control ########
  X <- try(as.matrix(X), silent = TRUE)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or coercible to a matrix.")
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", 1:ncol(X))
  if (!is.numeric(X)) stop("\"X\" only allows numerical or NA values.")
  if (anyNA(X)) {
    if (verbose) message("X contains NAs, dropping those observations.")
    mask_NA <- apply(X, 1, anyNA)
    X <- X[!mask_NA, , drop = FALSE]
  }
  if (nrow(X) <= ncol(X)) stop("Cannot deal with T <= N (after removing NAs), too few samples.")
  if (is.numeric(nu) && nu <= 2) stop("Non-valid value for nu (should be >2).")
  max_iter <- round(max_iter)
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)

  #
  # initialize parameters
  #
  alpha <- 1
  optimize_nu <- is.null(nu)
  if (optimize_nu) {
    nu    <- if (is.null(initial$nu)) 4 else initial$nu
  }
  if (nu == Inf) nu <- 1e5  # for numerical stability (for the Gaussian case)
  optimize_gamma <- is.null(gamma)
  if (optimize_gamma) {
    gamma_unscaled <- if (is.null(initial$gamma)) rep(0, N)  # sampleSkewness(X)
                      else initial$gamma
  } else
    gamma_unscaled <- gamma
  # # since mean = mu +  nu/(nu-2) * gamma:
  mu      <- if (is.null(initial$mu)) colMeans(X) - max(nu, 2.1)/(max(nu, 2.1)-2) * gamma_unscaled/alpha
             else initial$mu
  #mu      <- if (is.null(initial$mu)) colMeans(X) else initial$mu
  # since cov = nu/(nu-2) * scatter + 2*nu^2 / (nu-2)^2 / (nu-4) * gamma %o% gamma:
  scatter_unscaled <- if (is.null(initial$scatter)) (max(nu, 2.1)-2)/max(nu, 2.1) * cov(X)
                      else initial$scatter
  #scatter_unscaled <- if (is.null(initial$scatter)) cov(X) else initial$scatter

  # define tools for simplifying codes
  snapshot <- function() {
    if (ftol < Inf) list(mu = mu, gamma = gamma_unscaled/alpha, scatter = scatter_unscaled/alpha, nu = nu, alpha = alpha, obj = log_likelihood)
    else list(mu = mu, gamma = gamma_unscaled/alpha, scatter = scatter_unscaled/alpha, nu = nu, alpha = alpha)
  }


  #
  # loop
  #
  if (ftol < Inf)
    log_likelihood_record <- log_likelihood <- sum(dmvst(X = X, nu = nu, gamma = gamma_unscaled/alpha, mu = mu, scatter = scatter_unscaled/alpha))
  if (return_iterates) iterates_record <- list(snapshot())
  elapsed_times <- c(0)

  for (iter in 1:max_iter) {
    # record the current status
    start_time <- proc.time()[3]
    last_status <- snapshot()
    if (verbose) message(sprintf("iteration: %3d, objective: %.3f", iter, last_status$obj))

    # E-step ----------------------------------------
    # for the observed data
    expect <- Estep_mvst(X = X, nu = nu, gamma_unscaled = gamma_unscaled, mu = mu, scatter_unscaled = scatter_unscaled, alpha = alpha)

    # M-step ----------------------------------------
    # nu
    if (optimize_nu) {
      Q_nu <- function(nu) (nu/2)*sum(expect$E_logtau - log(alpha) - expect$E_tau/alpha) + T*((nu/2)*log(nu/2) - lgamma(nu/2))
      nu <- optimize(Q_nu, interval = c(getOption("nu_min"), getOption("nu_max")), maximum = TRUE)$maximum
    }

    # mu
    mu <- (colSums(X * expect$E_tau) - T*gamma_unscaled) / sum(expect$E_tau)

    # gamma
    if (optimize_gamma)
      gamma_unscaled <- (colSums(X) - T*mu) / sum(expect$E_invtau)

    # scatter
    X_ <- X - matrix(data = mu, nrow = T, ncol = N, byrow = TRUE)
    S <- - t(X_) %*% (X_ * expect$E_tau) / 2 +
      t(X_) %*% matrix(data = gamma_unscaled, nrow = T, ncol = N, byrow = TRUE) +
      - sum(expect$E_invtau) * gamma_unscaled %o% gamma_unscaled / 2
    scatter_unscaled <- - 2 * S / T
    scatter_unscaled <- (scatter_unscaled + t(scatter_unscaled)) / 2

    if (PXEM)
      alpha <- mean(expect$E_tau)


    # record the current the variables/loglikelihood if required
    if (ftol < Inf) {
      log_likelihood <- sum(dmvst(X = X, nu = nu, gamma = gamma_unscaled/alpha, mu = mu, scatter = scatter_unscaled/alpha))
      if (is.na(log_likelihood) || is.infinite(log_likelihood))
        stop("Error in computation of log-likelihood (nu = ", nu, ").")
      log_likelihood_old <- last_status$obj
      has_fun_converged <- abs(log_likelihood - log_likelihood_old) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_old))
    } else has_fun_converged <- TRUE
    if (ftol < Inf)
      log_likelihood_record <- c(log_likelihood_record, log_likelihood)
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    elapsed_times <- c(elapsed_times, proc.time()[3] - start_time)


    ## -------- stopping criterion --------
    have_params_converged <-
      all(abs(fnu(nu) - fnu(last_status$nu)) <= .5 * ptol * (abs(fnu(last_status$nu)) + abs(fnu(nu)))) &&
      all(abs(mu - last_status$mu)           <= .5 * ptol * (abs(mu) + abs(last_status$mu))) &&
      all(abs(gamma_unscaled/alpha - last_status$gamma)     <= .5 * ptol * (abs(gamma_unscaled/alpha) + abs(last_status$gamma))) &&
      all(abs(scatter_unscaled/alpha - last_status$scatter) <= .5 * ptol * (abs(scatter_unscaled/alpha) + abs(last_status$scatter)))
    if (is.na(have_params_converged) || is.na(has_fun_converged))
      browser()
    if (have_params_converged && has_fun_converged) break
  }


  # return results -------------
  gamma <- gamma_unscaled / alpha
  scatter <- scatter_unscaled / alpha
  # mean and cov matrix
  mean <- if (nu > 2)  mu +  nu/(nu-2) * gamma
          else NA
  cov <- if (nu > 4)  nu/(nu-2) * scatter + 2*nu^2 / (nu-2)^2 / (nu-4) * gamma %o% gamma
         else NA

  vars_to_be_returned <- list("mu"               = mu,
                              "gamma"            = gamma,
                              "scatter"          = scatter,
                              "nu"               = nu,
                              "mean"             = mean,
                              "cov"              = cov,
                              "converged"        = (iter < max_iter),
                              "num_iterations"   = iter,
                              "cpu_time"         = sum(elapsed_times))
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood_vs_iterations <- log_likelihood_record
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
    vars_to_be_returned$cpu_time_at_iter <- elapsed_times
  }
  return(vars_to_be_returned)
}



# sample estimator for skewness vector
sampleSkewness <- function(X) {
  sample_mean <- colMeans(X)
  X_ <- X - matrix(data = sample_mean, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  return(colMeans(X_^3) / (colMeans(X_^2)^1.5))
}


# log scaled probability density of give data if assumed to follow a GH MST distribution
#' @importFrom ghyp ghyp dghyp
dmvst_orig <- function(X, nu = 3, gamma = 1, mu = 0, scatter = 1) {
  model <- ghyp::ghyp(lambda = -nu/2, chi = nu, psi = 0, mu = mu, sigma = scatter, gamma = gamma)
  ghyp::dghyp(x = X, object = model, logvalue = TRUE)
}

dmvst <- function(X, nu = 3, gamma = 1, mu = 0, scatter = 1) {
  X <- as.matrix(X)
  N <- ncol(X)
  T <- nrow(X)
  Xc <- X - matrix(mu, T, N, byrow = TRUE)
  if (nu < 50000) {
    # model <- ghyp::ghyp(lambda = -nu/2, chi = nu, psi = 0, mu = mu, sigma = scatter, gamma = gamma)
    # return(ghyp::dghyp(x = X, object = model, logvalue = TRUE))
    scatter_inv <- solve(scatter)
    delta  <- as.numeric(sqrt(gamma %*% scatter_inv %*% gamma))
    if(delta == 0) {  # gamma = 0, student t case
      first_term <- lgamma((nu + N)/2) - lgamma(nu/2) -(N/2)*log(nu) - (N/2)*log(pi)-0.5* sum(log(eigen((scatter))$values) )
      second_term <- -((nu + N)/2) * log(1+(1/nu) * (rowSums(Xc * (Xc %*% scatter_inv))))
      return(first_term + second_term)
    } else {
      kappa <- sqrt(nu + rowSums(Xc * (Xc %*% scatter_inv)))
      first_term <-Xc %*% solve(scatter) %*% gamma - (N/2) * log(2*pi) - 0.5 * sum(log(eigen((scatter))$values) )
      second_term <- log(2) + (nu/2) * log(nu/2) - lgamma(nu/2)
      third_term <- -((nu+N)/2) * log(kappa/delta) + log_besselK(delta*kappa, -((nu + N)/2))
      return(first_term + second_term + third_term)
    }
  } else  # Gaussian case
    return(-0.5 * sum(log(eigen((2*pi*scatter))$values)) - 0.5 * rowSums(Xc * (Xc %*% solve(scatter))))
}



# expectation step of the EM algorithm for fitting a GH MST distribution
#' @importFrom numDeriv grad
Estep_mvst <- function(X, nu, gamma_unscaled, mu, scatter_unscaled, alpha) {
  N <- ncol(X)
  T <- nrow(X)

  gamma <- gamma_unscaled / alpha
  scatter <- scatter_unscaled / alpha
  scatter_inv <- solve(scatter)

  start_time <- proc.time()[3]  # record start time

  # three temporary values
  lambda <- (nu + N) / 2
  delta  <- as.numeric(sqrt(gamma %*% scatter_inv %*% gamma))
  Xc <- X - matrix(mu, T, N, byrow = TRUE)
  kappa <- sqrt(nu + rowSums(Xc * (Xc %*% scatter_inv)))
  # kappa  <- sqrt(nu + stats::mahalanobis(x = X, center = mu, cov = scatter_inv, inverted = TRUE))

  if (delta == 0) {
    E_tau <- (nu + N)/(kappa**2)
    E_invtau <- (kappa**2)/(nu + N -2)
    E_logtau <- digamma(lambda) - log((kappa**2)/2)
  } else {
    tmp <- besselK_ratio(delta * kappa, lambda)
    E_tau    <- (delta / kappa) * tmp
    #E_invtau  <- (kappa / delta) * 1/besselK_ratio(delta * kappa, lmd = lambda - 1)
    E_invtau <- (kappa / delta) * (tmp - (2*lambda)/delta/kappa)  # this saves computing bessel functions again

    if (lambda < 150) {
      dev_cal <- function(val) numDeriv::grad(func = function(lmd) log(besselK(x = val, nu = lmd, expon.scaled = TRUE)),
                                              x = lambda, method = "simple", method.args = list(eps = 1e-10))
      E_logtau <- log(delta / kappa) + sapply(delta * kappa, dev_cal)
    } else
      E_logtau <- log(delta / kappa) + log(besselK_ratio(delta * kappa, lambda))
  }

  # return
  list_to_return <- list("E_tau"    = E_tau * alpha,
                         "E_invtau" = E_invtau / alpha,
                         "E_logtau" = E_logtau + log(alpha))

  if (any(is.infinite(list_to_return$E_invtau)) ||
      any(is.nan(list_to_return$E_invtau)) ||
      any(is.nan(list_to_return$E_logtau))
      ) {
    message("Problem with the computation of E[tau], probably because of very small numbers in the evaluation of the bessel function.")
    browser()
  }

  return(list_to_return)
}



# https://www.researchgate.net/journal/Journal-of-Inequalities-and-Applications-1029-242X[…]mating-the-modified-Bessel-function-of-the-second-kind.pdf
besselK_ratio <- function(x, nu) {
  if (nu < 51)
    return(besselK(x = x, nu = nu + 1, expon.scaled = TRUE) / besselK(x = x, nu = nu, expon.scaled = TRUE))
  else {
    nu_i <- nu - floor(nu) + 50
    R_i <- besselK(x = x, nu = nu_i + 1, expon.scaled = TRUE) / besselK(x = x, nu = nu_i, expon.scaled = TRUE)
    while (nu_i != nu) {
      R_i <- 1/R_i + (2 * nu_i + 2)/x
      nu_i <- nu_i + 1
    }
    return(R_i)
  }
}


log_besselK <- function(x, nu) {
  if (nu <= 10 && nu >= -10) {
    return(log(besselK(x, nu)))
  } else if (nu >= 0){
    nu_i <-  nu - floor(nu) + 9
    K_nu_i <- besselK(x, nu = nu_i)    # K_{nu_i}(x)
    log_values <- log(K_nu_i)
    R_nu_i <- besselK(x, nu = nu_i + 1)/K_nu_i  # R_{nu_i}(x)
    log_values <- log_values + log(R_nu_i)
    nu_i <- nu_i + 1
    while (nu_i != nu) {
      R_nu_i <- 1/R_nu_i + 2 * nu_i /x
      log_values <- log_values + log(R_nu_i)
      nu_i <- nu_i + 1
    }
    return(log_values)
  } else {
    nu_i <- nu - floor(nu) - 9
    K_nu_i <- besselK(x, nu = nu_i)    # K_{nu_i}(x)
    log_values <- log(K_nu_i)
    R_nu_i <- K_nu_i/besselK(x, nu = nu_i - 1)  # R_{nu_i -1}(x)
    log_values <- log_values - log(R_nu_i)
    nu_i <- nu_i - 1
    while (nu_i != nu) {
      R_nu_i_inv <- R_nu_i - 2*nu_i/x
      R_nu_i <- 1/R_nu_i_inv   # R_{i -1}(x)
      log_values <- log_values - log(R_nu_i)
      nu_i <- nu_i - 1
    }
    return(log_values)
  }
}
