#' @title Estimate parameters of a multivariate (generalized hyperbolic) skewed t distribution to fit data
#'
#' @description Estimate parameters of a multivariate (generalized hyperbolic) skewed Student's t distribution to fit data,
#' namely, the location vector, the scatter matrix, the skewness vector, and the degrees of freedom.
#' The estimation is based on the maximum likelihood estimation (MLE) and the algorithm is
#' obtained from the expectation-maximization (EM) method.
#'
#' @details This function estimates the parameters of a (generalized hyperbolic) multivariate Student's t distribution (\code{mu},
#'          \code{scatter}, \code{gamme} and \code{nu}) to fit the data via the expectation-maximization (EM) algorithm.
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{nu}: default is \code{4},}
#'                         \item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{gamme}: default is the sample skewness vector,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         }
#' @param max_iter Integer indicating the maximum number of iterations for the iterative estimation
#'                 method (default is \code{100}).
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
#' @return A list containing possibly the following elements:
#'         \item{\code{mu}}{Location vector estimate.}
#'         \item{\code{gamma}}{Skewness vector estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{cov}}{Covarance matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has been reached (\code{FALSE}).}
#'         \item{\code{num_iterations}}{Number of iterations executed.}
#'         \item{\code{cpu_time}}{Elapsed overall CPU time.}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu}, \code{scatter}, \code{nu},
#'                                       and possibly \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (if \code{return_iterates = TRUE}).}
#'         \item{\code{cpu_time_at_iter}}{Elapsed CPU time at each iteration (if \code{return_iterates = TRUE}).}
#'
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @seealso \code{\link{fit_mvt}}
#'
#' @references
#' Aas, Kjersti and Ingrid Hobæk Haff. "The generalized hyperbolic skew student’st-distribution,"
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
#' # fit GH Skew t model
#' fit_mvst(X)
#'
#' @importFrom stats optimize
#' @export
fit_mvst <- function(X, initial = NULL, max_iter = 100, ptol = 1e-3, ftol = Inf,
                     PXEM = TRUE, return_iterates = FALSE, verbose = FALSE) {

  T <- nrow(X)
  N <- ncol(X)

  # initialize parameters
  alpha   <- 1
  nu      <- if (is.null(initial$nu))      4                   else initial$nu
  mu      <- if (is.null(initial$mu))      colMeans(X)         else initial$mu
  gamma   <- if (is.null(initial$gamma))   sampleSkewness(X)   else initial$gamma
  scatter <- if (is.null(initial$scatter)) cov(X)              else initial$scatter

  # define tools for simplifying codes
  snapshot <- function() {
    snap <- list(nu = nu, mu = mu, gamma = gamma/alpha, scatter = scatter/alpha, alpha = alpha)
    if (verbose || ftol < Inf) snap$obj = sum(dST(X = X, nu = nu, gamma = gamma/alpha, mu = mu, scatter = scatter/alpha))
    return(snap)
  }

  # let's loop
  if (return_iterates) iterates_record <- list(snapshot())
  elapsed_times <- c(0)

  for (iter in 1:max_iter) {

    start_time <- proc.time()[3]  # record start time

    last_status <- snapshot()

    if (verbose) message(sprintf("iteration: %3d, objective: %.3f", iter, last_status$obj))

    # E-step ----------------------------------------
    # for the observed data
    # browser()
    expect <- Estep_mst(X = X, nu = nu, gamma = gamma, mu = mu, scatter = scatter, alpha = alpha)

    # M-step ----------------------------------------
    # nu
    Q_nu <- function(nu) (nu/2)*sum(expect$E_logtau - log(alpha) - expect$E_tau/alpha) + T*((nu/2)*log(nu/2) - log(base::gamma(nu/2)))
    nu <- optimize(Q_nu, interval = c(2 + 1e-12, 50), maximum = TRUE)$maximum

    # mu
    mu <- (colSums(X * expect$E_tau) - T*gamma) / sum(expect$E_tau)

    # gamma
    gamma <- (colSums(X) - T*mu) / sum(expect$E_invtau)

    # scatter
    X_ <- X - matrix(data = mu, nrow = T, ncol = N, byrow = TRUE)
    S <- - t(X_) %*% (X_ * expect$E_tau) / 2 +
      t(X_) %*% matrix(data = gamma, nrow = T, ncol = N, byrow = TRUE) +
      - sum(expect$E_invtau) * cbind(gamma) %*% rbind(gamma) / 2
    scatter <- - 2 * S / T
    scatter <- (scatter + t(scatter)) / 2

    if (PXEM) alpha <- mean(expect$E_tau)

    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()

    have_params_converged <-
      all(abs(fnu(nu) - fnu(last_status$nu)) <= .5 * ptol * (abs(fnu(last_status$nu)) + abs(fnu(nu)))) &&
      all(abs(mu - last_status$mu)           <= .5 * ptol * (abs(mu) + abs(last_status$mu))) &&
      all(abs(gamma - last_status$gamma)     <= .5 * ptol * (abs(gamma) + abs(last_status$gamma))) &&
      all(abs(scatter - last_status$scatter) <= .5 * ptol * (abs(scatter) + abs(last_status$scatter)))

    if (ftol < Inf) {
      log_likelihood_current <- snapshot()$obj
      log_likelihood_old     <- last_status$obj
      has_fun_converged <- abs(log_likelihood_current - log_likelihood_old) <= .5 * ftol * (abs(log_likelihood_current) + abs(log_likelihood_old))
    } else has_fun_converged <- TRUE

    elapsed_times <- c(elapsed_times, proc.time()[3] - start_time)

    if (have_params_converged && has_fun_converged) break
  }


  # return results -------------
  gamma <- gamma / alpha
  scatter <- scatter / alpha
  # cov matrix
  if (nu > 4)
    Sigma_cov <- nu/(nu-2) * scatter + 2*nu^2 / (nu-2)^2 / (nu-4) * cbind(gamma) %*% rbind(gamma)
  else
    Sigma_cov <- NA

  vars_to_be_returned <- list("mu"               = mu,
                              "gamma"            = gamma,
                              "scatter"          = scatter,
                              "cov"              = Sigma_cov,
                              "nu"               = nu,
                              "converged"        = (iter < max_iter),
                              "num_iterations"   = iter,
                              "cpu_time"         = sum(elapsed_times))
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
    vars_to_be_returned$cpu_time_at_iterelapsed_times
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
dST <- function(X, nu = 3, gamma = 1, mu = 0, scatter = 1) {
  model <- ghyp::ghyp(lambda = -nu/2, chi = nu, psi = 0, mu = mu, sigma = scatter, gamma = gamma)
  ghyp::dghyp(x = X, object = model, logvalue = TRUE)
}


# expectation step of the EM algorithm for fitting a GH MST distribution
#' @importFrom numDeriv grad
#' @importFrom stats mahalanobis
Estep_mst <- function(X, nu, gamma, mu, scatter, alpha) {
  N <- ncol(X)
  T <- nrow(X)

  gamma <- gamma / alpha
  scatter <- scatter / alpha

  scatter_inv <- solve(scatter)

  start_time <- proc.time()[3]  # record start time

  # three temporary values
  lambda <- (nu + N) / 2
  delta  <- as.numeric(sqrt(gamma %*% scatter_inv %*% gamma))
  kappa  <- sqrt(nu + mahalanobis(x = X, center = mu, cov = scatter_inv, inverted = TRUE))

  E_tau    <- (delta / kappa) * besselK(x = delta * kappa, nu = lambda + 1, expon.scaled = TRUE) / besselK(x = delta * kappa, nu = lambda, expon.scaled = TRUE)
  E_invtau <- (kappa / delta) * besselK(x = delta * kappa, nu = lambda - 1, expon.scaled = TRUE) / besselK(x = delta * kappa, nu = lambda, expon.scaled = TRUE)

  dev_cal <- function(val) numDeriv::grad(func = function(lmd) log(besselK(x = val, nu = lmd, expon.scaled = TRUE)), x = lambda, method = "simple", method.args = list(eps = 1e-10))
  E_logtau <- log(delta / kappa) + sapply(delta*kappa, dev_cal)

  list("E_tau"    = E_tau * alpha,
       "E_invtau" = E_invtau / alpha,
       "E_logtau" = E_logtau + log(alpha))
}

