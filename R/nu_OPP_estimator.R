#' @title Estimate the degrees of freedom of a heavy-tailed t distribution based on the OPP estimator
#'
#' @description This function estimates the degrees of freedom of a heavy-tailed \eqn{t} distribution based on
#'              the OPP estimator from paper [Ollila-Palomar-Pascal, TSP2021, Alg. 1].
#'              Traditional nonparametric methods or likelihood methods provide erratic estimations of
#'              the degrees of freedom unless the number of observations is very large.
#'              The POP estimator provides a stable estimator based on random matrix theory.
#'              A number of different versions are provided, but the default POP method will most likely
#'              be the desired choice.
#'
#' @param var_X Vector with the sample variance of the columns of the data matrix.
#' @param trace_scatter Trace of the scatter matrix.
#' @param r2 Vector containing the values of \code{diag( Xc \%*\% inv(scatter) \%*\% t(Xc) )},
#'           where \code{Xc} is the centered data matrix.
#' @param method String indicating the version of the OPP estimator (default is just \code{"OPP"}).
#'               Other option is the variation: \code{"OPP-harmonic"}.
#'
#' @return Estimated value of the degrees of freedom \code{nu} of a heavy-tailed \eqn{t} distribution.
#'
#' @author Esa Ollila, Frédéric Pascal, and Daniel P. Palomar
#'
#' @references
#'
#' Esa Ollila, Daniel P. Palomar, and Frédéric Pascal, "Shrinking the Eigenvalues of M-estimators of Covariance Matrix,"
#' IEEE Trans. on Signal Processing, vol. 69, pp. 256-269, Jan. 2021. <https://doi.org/10.1109/TSP.2020.3043952>
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
#' nu_OPP_estimator(var_X = 1/(T-1)*colSums(Xc^2), trace_scatter = sum(diag(Sigma_true)))
#'
#' # usage #2
#' r2 <- rowSums(Xc * (Xc %*% solve(Sigma_true)))
#' nu_OPP_estimator(var_X = 1/(T-1)*colSums(Xc^2), trace_scatter = sum(diag(Sigma_true)), method = "OPP-harmonic", r2 = r2)
#'
#' @export
nu_OPP_estimator <- function(var_X, trace_scatter, r2,
                             method = c("OPP", "OPP-harmonic")) {
  N <- length(var_X)

  eta <- switch(match.arg(method),
                "OPP"          = sum(var_X)/trace_scatter,
                "OPP-harmonic" = sum(var_X)/trace_scatter * mean(N/r2),
                stop("OPP method unknown.")
                )
  nu <- 2*eta/(eta - 1)

  # cap nu and return
  nu <- min(getOption("nu_max"), max(getOption("nu_min"), nu))
  return(nu)
}

