#' @title Estimate the degrees of freedom of a heavy-tailed t distribution based on the POP estimator
#'
#' @description This function estimates the degrees of freedom of a heavy-tailed \eqn{t} distribution based on
#'              the POP estimator from paper [Pascal-Ollila-Palomar, EUSIPCO2021, Alg. 1].
#'              Traditional nonparametric methods or likelihood methods provide erratic estimations of
#'              the degrees of freedom unless the number of observations is very large.
#'              The POP estimator provides a stable estimator based on random matrix theory.
#'              A number of different versions are provided, but the default POP method will most likely
#'              be the desired choice.
#'
#' @param Xc Centered data matrix (with zero mean) containing the multivariate time series (each column is one time series).
#' @param N Number of variables (columns of data matrix) in the multivariate time series.
#' @param T Number of observations (rows of data matrix) in the multivariate time series.
#' @param nu Current estimate of the degrees of freedom of the \eqn{t} distribution.
#' @param Sigma Current estimate of the scatter matrix.
#' @param alpha Value for the acceleration technique (cf. \code{\link{fit_mvt}()}).
#' @param r2 Vector containing the values of \code{diag( Xc \%*\% inv(scatter) \%*\% t(Xc) )}.
#' @param method String indicating the version of the POP estimator (default is just \code{"POP"} and should work well in all cases).
#'               Other versions include: \code{"POP-approx-1"}, \code{"POP-approx-2"}, \code{"POP-approx-3"},
#'               \code{"POP-approx-4"}, \code{"POP-exact"}, \code{"POP-sigma-corrected"}, \code{"POP-sigma-corrected-true"}.
#'
#' @return Estimated value of the degrees of freedom \code{nu} of a heavy-tailed \eqn{t} distribution.
#'
#' @author Frédéric Pascal, Esa Ollila, and Daniel P. Palomar
#'
#' @references
#'
#' Frédéric Pascal, Esa Ollila, and Daniel P. Palomar, "Improved estimation of the degree of freedom parameter of
#' multivariate t-distribution," in Proc. European Signal Processing Conference (EUSIPCO), Dublin, Ireland, Aug. 23-27, 2021.
#' <https://doi.org/10.23919/EUSIPCO54536.2021.9616162>
#'
#' @examples
#' library(mvtnorm)       # to generate heavy-tailed data
#' library(fitHeavyTail)
#'
#' # parameters
#' N <- 5
#' T <- 100
#' nu_true <- 4           # degrees of freedom
#' mu_true <- rep(0, N)   # mean vector
#' Sigma_true <- diag(N)  # scatter matrix
#'
#' # generate data
#' X <- rmvt(n = T, sigma = Sigma_true, delta = mu_true, df = nu_true)  # generate Student's t data
#' mu <- colMeans(X)
#' Xc <- X - matrix(mu, T, N, byrow = TRUE)    # center data
#'
#' # usage #1
#' nu_POP_estimator(Xc = Xc, nu = 10, Sigma = Sigma_true)
#'
#' # usage #2
#' r2 <- rowSums(Xc * (Xc %*% solve(Sigma_true)))
#' nu_POP_estimator(r2 = r2, nu = 10, N = N)
#'
#' # usage #3
#' nu_POP_estimator(r2 = r2, nu = 10, N = N, method = "POP-approx-1")
#'
#' @importFrom stats uniroot
#' @importFrom mvtnorm rmvt
#' @export
nu_POP_estimator <- function(Xc = NULL, N = NULL, T = NULL,
                             Sigma = NULL, nu = NULL,
                             r2 = NULL,
                             method = c("POP",
                                        "POP-approx-1", "POP-approx-2", "POP-approx-3", "POP-approx-4", "POP-exact",
                                        "POP-sigma-corrected", "POP-sigma-corrected-true"),
                             alpha = 1) {
  # methods
  method <- match.arg(method)
  if (method == "POP")
    method <- "POP-approx-2"  # default POP method

  # get N and T automatically from other variables if possible
  if (is.null(T)) {
    if (!is.null(r2))
      T <- length(r2)
    else if (!is.null(Xc))
      T <- nrow(Xc)
  }
  if (is.null(N)) {
    if (!is.null(Xc))
      N <- ncol(Xc)
    else if (!is.null(Sigma))
      N <- ncol(Sigma)
  }
  if (is.null(T) || is.null(N))
    stop("POP estimator needs to know both N and T.")
  if (is.null(nu))
    stop("POP estimator needs the previous estimation of nu.")
  if (is.null(r2)) {
    if (is.null(Sigma))
      stop("POP estimator needs Sigma to compute r2.")
    r2 <- rowSums(Xc * (Xc %*% solve(Sigma)))
  }

  # different versions of POP estimator
  theta <- switch(method,
                  "POP-approx-1" = {
                    u <- (N + nu)/(nu + r2)
                    r2i <- r2/(1 - r2*u/T)
                    (1 - N/T) * sum(r2i)/T/N
                  },
                  "POP-approx-2" = {  # <<<--------    default POP method
                    u <- (N + nu)/(nu + r2*T/(T-1))
                    r2i <- r2/(1 - r2*u/T)
                    (1 - N/T) * sum(r2i)/T/N
                  },
                  "POP-approx-3" = {  # not good
                    Sigma_ <- (1/T) * crossprod(Xc)
                    r2_ <- rowSums(Xc * (Xc %*% solve(Sigma_)))
                    r2i <- r2_/(1 - r2_/T)
                    (1 - N/T) * sum(r2i)/T/N
                  },
                  "POP-approx-4" = {  # not good
                    r2i <- r2/(1 - r2/T)
                    (1 - N/T) * sum(r2i)/T/N
                  },
                  "POP-exact" = {
                    u <- (N + nu)/(nu + r2)
                    r2i <- vector("numeric", length = T)
                    for (ii in 1:T) {
                      # hardcoded iteration #1
                      ui <- u
                      Sigmai <- (1/T/alpha) * crossprod(sqrt(ui[-ii]) * Xc[-ii, ])
                      r2i_ <- rowSums(Xc * (Xc %*% solve(Sigmai)))
                      # hardcoded iteration #2
                      ui <- (N + nu) / (nu + r2i_)
                      Sigmai <- (1/T/alpha) * crossprod(sqrt(ui[-ii]) * Xc[-ii, ])
                      r2i_ <- rowSums(Xc * (Xc %*% solve(Sigmai)))
                      # hardcoded iteration #3
                      ui <- (N + nu) / (nu + r2i_)
                      Sigmai <- (1/T/alpha) * crossprod(sqrt(ui[-ii]) * Xc[-ii, ])
                      r2i[ii] <- as.numeric(Xc[ii, ] %*% solve(Sigmai, Xc[ii, ]))
                    }
                    (1 - N/T) * sum(r2i)/T/N
                  },
                  "POP-sigma-corrected" = {
                    u <- (N + nu)/(nu + r2)
                    r2i <- r2/(1 - r2*u/T)
                    theta <- (1 - N/T) * sum(r2i)/T/N
                    # correction with sigma
                    psi <- function(t) (N + nu)/(nu + t) * t
                    F <- function(sigma) mean(psi(r2/sigma)) - N  # should I use r2i in here?
                    sigma <- uniroot(F, lower = 0.1, upper = 1000)$root
                    theta <- theta * sigma
                    theta
                  },
                  "POP-sigma-corrected-true" = {
                    u <- (N + nu)/(nu + r2)
                    r2i <- r2/(1 - r2*u/T)
                    theta <- (1 - N/T) * sum(r2i)/T/N
                    # correction with sigma
                    T_ <- 10000  #10000
                    X_ <- rmvt(n = T_, delta = rep(0, N), sigma = diag(N), df = nu)  # heavy-tailed data
                    r2_ <- rowSums(X_^2)
                    u_ <- (N + nu)/(nu + r2_)
                    r2i_ <- r2_/(1 - r2_*u_/T_)
                    psi <- function(t) (N + nu)/(nu + t) * t
                    F <- function(sigma) mean(psi(r2_/sigma)) - N
                    sigma <- uniroot(F, lower = 0.1, upper = 1000)$root
                    theta <- theta * sigma
                    theta
                  },
                  "theta-1a" = sum(r2)/T/N,
                  "theta-1b" = {
                    r2 <- rowSums(Xc * (Xc %*% solve(Sigma)))
                    T/(T-1)*sum(r2)/T/N
                  },
                  "theta-2b" = {
                    u <- (N + nu)/(nu + r2)
                    r2i <- r2/(1 - r2*u/T)
                    (T - N + 2)*T/(T - 1)/(T + 1) * sum(r2i)/T/N  #(1 - (N + 2)/T) * sum(r2i)/T/N
                  },
                  stop("POP method unknown.")
                )
  nu <- 2*theta/(theta - 1)

  # cap nu and return
  nu <- min(getOption("nu_max"), max(getOption("nu_min"), nu))
  return(nu)
}

