#
# lower and upper bounds for the optimization of nu
#
.nu_min <- 2.5
.nu_max <- 100


#' @title Estimate parameters of a multivariate Student's t distribution to fit data
#'
#' @description Estimate parameters of a multivariate Student's t distribution to fit data,
#' namely, the mean vector, the covariance matrix, the scatter matrix, and the degrees of freedom.
#' The data can contain missing values denoted by NAs.
#' It can also consider a factor model structure on the covariance matrix.
#' The estimation is based on the maximum likelihood estimation (MLE) and the algorithm is
#' obtained from the expectation-maximization (EM) method.
#'
#' @details This function estimates the parameters of a multivariate Student's t distribution (\code{mu},
#'          \code{cov}, \code{scatter}, and \code{nu}) to fit the data via the expectation–maximization (EM) algorithm.
#'          The data matrix \code{X} can contain missing values denoted by NAs.
#'          The estimation of \code{nu} if very flexible: it can be directly passed as an argument (without being estimated),
#'          it can be estimated with several one-shot methods (namely, \code{"kurtosis"}, \code{"MLE-diag"},
#'          \code{"MLE-diag-resampled"}), and it can also be iteratively estimated with the other parameters via the EM
#'          algorithm.
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param na_rm Logical value indicating whether to remove observations with some NAs (default) or not, in which
#'              case they will be imputed at a higher computational cost.
#' @param nu Degrees of freedom of the \eqn{t} distribution. Either a number (\code{>2}) or a string indicating the
#'           method to compute it:
#'           \itemize{\item{\code{"kurtosis"}: based on the kurtosis obtained from the sampled moments;}
#'                    \item{\code{"MLE-diag"}: based on the MLE assuming a diagonal sample covariance;}
#'                    \item{\code{"MLE-diag-resampled"}: method "MLE-diag" resampled for better stability;}
#'                    \item{\code{"iterative"}: iterative estimation with the rest of the parameters via the EM algorithm.}}
#' @param nu_iterative_method String indicating the method for iteratively estimating \code{nu} (in case \code{nu = "iterative"}):
#'                  \itemize{\item{\code{"ECM"}: maximization of the Q function;}
#'                           \item{\code{"ECME"}: maximization of the log-likelihood function;}
#'                           \item{\code{"ECME-diag"}: maximization of the log-likelihood function assuming
#'                                                     a digonal scatter matrix (default method).}}
#'                  This argument is used only when there are no NAs in the data and no factor model is chosen.
#' @param initial List of initial values of the parameters for the iterative EM estimation method (in case \code{nu = "iterative"}).
#'                Possible elements include:
#'                \itemize{\item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{cov}: default is the data sample covariance matrix,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         \item{\code{nu}: can take the same values as argument \code{nu}, default is \code{4},}
#'                         \item{\code{B}: default is the top eigenvectors of \code{initial$cov}
#'                                                   multiplied by the sqrt of the eigenvalues,}
#'                         \item{\code{psi}: default is
#'                                          \code{diag(initial$cov - initial$B \%*\% t(initial$B)).}}}
#' @param factors Integer indicating number of factors (default is \code{ncol(X)}, so no factor model assumption).
#' @param max_iter Integer indicating the maximum number of iterations for the iterative estimation
#'                 method (default is \code{100}).
#' @param ptol Positive number indicating the relative tolerance for the change of the variables
#'             to determine convergence of the iterative method (default is \code{1e-3}).
#' @param ftol Positive number indicating the relative tolerance for the change of the log-likelihood
#'             value to determine convergence of the iterative method (default is \code{Inf}, so it is
#'             not active). Note that using this argument might have a computational cost as a convergence
#'             criterion due to the computation of the log-likelihood (especially when \code{X} is high-dimensional).
#' @param return_iterates Logical value indicating whether to record the values of the parameters (and possibly the
#'                        log-likelihood if \code{ftol < Inf}) at each iteration (default is \code{FALSE}).
#' @param verbose Logical value indicating whether to allow the function to print messages (default is \code{FALSE}).
#'
#' @return A list containing possibly the following elements:
#'         \item{\code{mu}}{Mean vector estimate.}
#'         \item{\code{cov}}{Covariance matrix estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has been reached (\code{FALSE}).}
#'         \item{\code{num_iterations}}{Number of iterations executed.}
#'         \item{\code{cpu_time}}{Elapsed CPU time.}
#'         \item{\code{B}}{Factor model loading matrix estimate according to \code{cov = (B \%*\% t(B) + diag(psi)}
#'                         (only if factor model requested).}
#'         \item{\code{psi}}{Factor model idiosynchratic variances estimates according to \code{cov = (B \%*\% t(B) + diag(psi)}
#'                           (only if factor model requested).}
#'         \item{\code{log_likelihood}}{Value of log-likelihood after converge of the estimation algorithm (if \code{ftol < Inf}).}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu}, \code{scatter}, \code{nu},
#'                                       and possibly \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (if \code{return_iterates = TRUE}).}
#'
#' @author Daniel P. Palomar and Rui Zhou
#'
#' @seealso \code{\link{fit_Tyler}} and \code{\link{fit_Cauchy}}
#'
#' @references
#' Chuanhai Liu and Donald B. Rubin, “ML estimation of the t-distribution using EM and its extensions, ECM and ECME,”
#' Statistica Sinica (5), pp. 19-39, 1995.
#'
#' Rui Zhou, Junyan Liu, Sandeep Kumar, and Daniel P. Palomar, "Robust factor analysis parameter estimation,"
#' Lecture Notes in Computer Science (LNCS), 2019. <https://arxiv.org/abs/1909.12530>
#'
#' @examples
#' library(mvtnorm)       # to generate heavy-tailed data
#' library(fitHeavyTail)
#'
#' X <- rmvt(n = 1000, df = 6)  # generate Student's t data
#' fit_mvt(X)
#'
#' @importFrom stats optimize
#' @export
fit_mvt <- function(X, na_rm = TRUE,
                    nu = c("kurtosis", "MLE-diag", "MLE-diag-resampled", "iterative"),
                    nu_iterative_method = c("ECME-diag", "ECME", "ECM", "ECME-cov"),
                    initial = NULL, factors = ncol(X),
                    max_iter = 100, ptol = 1e-3, ftol = Inf,
                    return_iterates = FALSE, verbose = FALSE) {
  ####### error control ########
  X <- try(as.matrix(X), silent = TRUE)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or coercible to a matrix.")
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", 1:ncol(X))
  if (!is.numeric(X)) stop("\"X\" only allows numerical or NA values.")
  if (anyNA(X))
    if (na_rm || ncol(X) == 1) {
      if (verbose) message("X contains NAs, dropping those observations.")
      mask_NA <- apply(X, 1, anyNA)
      X <- X[!mask_NA, , drop = FALSE]
    }
  if (nrow(X) <= ncol(X)) stop("Cannot deal with T <= N (after removing NAs), too few samples.")
  if (is.numeric(nu) && nu <= 2) stop("Non-valid value for nu.")
  factors <- round(factors)
  if (factors < 1 || factors > ncol(X)) stop("\"factors\" must be no less than 1 and no more than column number of \"X\".")
  max_iter <- round(max_iter)
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)
  X_has_NA <- anyNA(X)
  FA_struct <- (factors != N)
  if (!is.numeric(nu)) nu <- match.arg(nu)
  nu_iterative_method <- match.arg(nu_iterative_method)
  optimize_nu <- (nu == "iterative")

  # initialize all parameters
  start_time <- proc.time()[3]
  if (optimize_nu) {  # initial point
    nu <- if (!is.null(initial$nu)) initial$nu
          else nu <- 4  # default initial point
  }
  if (!is.numeric(nu)) {
    nu <- switch(nu,
                 "kurtosis"           = nu_from_kurtosis(X),
                 "MLE-diag"           = nu_mle(X, method = "MLE-mv-diagcov"),
                 "MLE-diag-resampled" = nu_mle(X, method = "MLE-mv-diagcov-resampled"),
                 stop("Method to estimate nu unknown."))
    if (verbose) message(sprintf("Automatically setting nu = %.2f", nu))
  } else if (nu == Inf) nu <- 1e15  # for numerical stability (for the Gaussian case)
  mu <- if (is.null(initial$mu)) colMeans(X, na.rm = TRUE) else initial$mu
  Sigma <- if (is.null(initial$cov)) (nu-2)/nu * var(X, na.rm = TRUE) else (nu-2)/nu * initial$cov
  if (!is.null(initial$scatter)) Sigma <- initial$scatter
  if (FA_struct) {  # Sigma is the scatter matrix, not the covariance matrix
    Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
    B <- if (is.null(initial$B)) Sigma_eigen$vectors[, 1:factors] %*% diag(sqrt(Sigma_eigen$values[1:factors]), factors)
         else initial$B
    psi <- if (is.null(initial$psi)) pmax(0, diag(Sigma) - diag(B %*% t(B)))
           else initial$psi
    Sigma <- B %*% t(B) + diag(psi, N)
  }
  if (ftol < Inf) log_likelihood <- ifelse(X_has_NA,
                                           dmvt_withNA(X = X, delta = mu, sigma = Sigma / alpha, df = nu),
                                           sum(mvtnorm::dmvt(X, delta = mu, sigma = Sigma, df = nu, log = TRUE, type = "shifted")))
  alpha <- 1  # an extra variable for PX-EM acceleration

  # aux function to save iterates
  snapshot <- function() {
    if (ftol < Inf) list(mu = mu, scatter = Sigma, nu = nu, log_likelihood = log_likelihood)
    else list(mu = mu, scatter = Sigma, nu = nu)
  }

  # loop
  if (return_iterates) iterates_record <- list(snapshot())
  for (iter in 1:max_iter) {
    # record the current status
    mu_old    <- mu
    Sigma_old <- Sigma
    nu_old    <- nu
    if (ftol < Inf) log_likelihood_old <- log_likelihood

    ## -------------- E-step --------------
    if (X_has_NA || FA_struct)
      Q <- Estep(mu, Sigma, psi, nu, X)
    else {
      X_ <- X - matrix(mu, T, N, byrow = TRUE)
      delta2 <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
      E_tau <- (nu + N) / (nu + delta2)
      ave_E_tau <- mean(E_tau)
      ave_E_tau_X <- (1/T)*as.vector(E_tau %*% X)
    }

    ## -------------- M-step --------------
    # update mu, alpha, nu
    if (X_has_NA || FA_struct) {
      mu <- Q$ave_E_tau_X / Q$ave_E_tau
      alpha <- Q$ave_E_tau
      S <- Q$ave_E_tau_XX - cbind(mu) %*% rbind(Q$ave_E_tau_X) - cbind(Q$ave_E_tau_X) %*% rbind(mu) + Q$ave_E_tau * cbind(mu) %*% rbind(mu)
      S <- S / alpha
      if (FA_struct) {
        B   <- optB(S = S, factors = factors, psi_vec = psi)
        psi <- pmax(0, diag(S - B %*% t(B)))
        Sigma <- B %*% t(B) + diag(psi, N)
      } else
        Sigma <- S
      if (optimize_nu) {
        Q_nu <- function(nu) - (nu/2)*log(nu/2) + lgamma(nu/2) - (nu/2)*(Q$ave_E_logtau - Q$ave_E_tau)
        nu <- optimize(Q_nu, interval = c(2 + 1e-12, 100))$minimum
      }
    } else {
      mu <- ave_E_tau_X / ave_E_tau
      alpha <- ave_E_tau  # acceleration
      X_ <- X - matrix(mu, T, N, byrow = TRUE)  # this is slower: sweep(X, 2, FUN = "-", STATS = mu)  #X_ <- X - rep(mu, each = TRUE)  # this is wrong?
      ave_E_tau_XX <- (1/T) * crossprod(sqrt(E_tau) * X_)  # (1/T) * t(X_) %*% diag(E_tau) %*% X_
      Sigma <- ave_E_tau_XX / alpha
      if (optimize_nu)
        nu <- switch(nu_iterative_method,
                     "ECM" = {  # based on minus the Q function of nu
                       #nu_regcoef = 1e2
                       #nu_target <- nu_from_kurtosis(X)
                       ave_E_log_tau_minus_E_tau <- digamma((N+nu)/2) - log((N+nu)/2) + mean(log(E_tau) - E_tau)  # equal to Q$ave_E_logtau - Q$ave_E_tau
                       Q_nu <- function(nu) - (nu/2)*log(nu/2) + lgamma(nu/2) - (nu/2)*ave_E_log_tau_minus_E_tau
                                              #+ nu_regcoef * (1/T)*(nu/(nu-2) - nu_target/(nu_target-2))^2
                       optimize(Q_nu, interval = c(.nu_min, .nu_max))$minimum
                       },
                     "ECME" = {  # based on minus log-likelihood of nu with mu and sigma fixed to mu[k+1] and sigma[k+1]
                       #nu_mle(Xc = X_, method = "MLE-mv-scat", Sigma_scatter = Sigma)
                       delta2 <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
                       negLL <- function(nu) ((N + nu)/2)*sum(log(nu + delta2)) - T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                       optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
                       },
                     "ECME-diag" = {  # based on minus log-likelihood of nu with mu and sigma fixed to mu[k+1] and sigma[k+1]
                       #nu_mle(Xc = X_, method = "MLE-mv-diagscat", Sigma_scatter = Sigma)
                       delta2 <- rowSums(X_^2 / matrix(diag(Sigma), T, N, byrow = TRUE))  # this amounts to using a diagonal cov matrix
                       negLL <- function(nu) ((N + nu)/2)*sum(log(nu + delta2)) - T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                       optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
                     },
                     "ECME-cov" = {  # this variation is worse than ECME
                       Sigma_cov <- nu/(nu-2)*Sigma
                       #nu_mle(Xc = X_, method = "MLE-mv-cov", Sigma_cov = Sigma_cov)  # this is without nu_target
                       delta2_cov <- rowSums(X_ * (X_ %*% inv(Sigma_cov)))
                       negLL <- function(nu) (N*T)/2*log((nu-2)/nu) + ((N + nu)/2)*sum(log(nu + nu/(nu-2)*delta2_cov)) -
                         T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                       optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
                     },
                     stop("Method to estimate nu unknown."))
    }

    ## -------- stopping criterion --------
    have_params_converged <-
      all(abs(mu - mu_old)       <= .5 * ptol * (abs(mu_old) + abs(mu))) &&
      abs(fnu(nu) - fnu(nu_old)) <= .5 * ptol * (abs(fnu(nu_old)) + abs(fnu(nu))) &&
      all(abs(Sigma - Sigma_old) <= .5 * ptol * (abs(Sigma_old) + abs(Sigma)))

    if (ftol < Inf) {
      log_likelihood  <- dmvt_withNA(X = X, delta = mu, sigma = Sigma, df = nu)
      has_fun_converged <- abs(log_likelihood - log_likelihood_old) <= .5 * ftol * (abs(log_likelihood) + abs(log_likelihood_old))
    } else has_fun_converged <- TRUE
    # record the current the variables/loglikelihood if required
    names(mu) <- colnames(Sigma) <- rownames(Sigma) <- colnames(X)  # add names first
    if (return_iterates) iterates_record[[iter + 1]] <- snapshot()
    if (have_params_converged && has_fun_converged) break
  }
  elapsed_time <- proc.time()[3] - start_time
  if (verbose) message(sprintf("Number of iterations for mvt estimation = %d\n", iter))

  ## -------- return variables --------
  #Sigma <- T/(T-1) * Sigma  # unbiased estimator
  vars_to_be_returned <- list("mu"             = mu,
                              "cov"            = if (nu > 2) nu/(nu-2) * Sigma else NA,
                              "scatter"        = Sigma,
                              "nu"             = nu,
                              "converged"      = (iter < max_iter),
                              "num_iterations" = iter,
                              "cpu_time"       = elapsed_time)
  # if (!optimize_nu) {
  #   kappa <- scaling_fitting_ka_with_b(a = diag(Sigma), b = apply(X^2, 2, mean, trim = max(1/T, 0.03)))
  #   vars_to_be_returned$cov <- kappa * Sigma
  # }
  if (FA_struct) {
    rownames(B) <- names(psi) <- colnames(X)
    colnames(B) <- paste0("factor-", 1:ncol(B))
    vars_to_be_returned$B   <- sqrt(nu/(nu-2)) * B
    vars_to_be_returned$psi <- nu/(nu-2) * psi
  }
  if (ftol < Inf)
    vars_to_be_returned$log_likelihood <- log_likelihood
  if (return_iterates) {
    names(iterates_record) <- paste("iter", 0:(length(iterates_record)-1))
    vars_to_be_returned$iterates_record <- iterates_record
  }
  return(vars_to_be_returned)
}




##
## -------- Auxiliary functions --------
##


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
# optimize_nu <- function(y, tol = 1e-4) {
#   if (y < 0) stop("y must be positive!")
#   L <- 1e-5
#   U <- 1000
#
#   while ((log(U)-digamma(U)) > y) U <- U*2
#   while ((log(L)-digamma(L)) < y) L <- L*2
#
#   while (1) {
#     mid <- (L + U) / 2
#     tmp <- log(mid)-digamma(mid) - y
#
#     if (abs(tmp) < tol) break
#     if (tmp > 0) L <- mid else U <- mid
#   }
#
#   return(2*mid)
# }


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


# estimate nu via the kurtosis of each variable
excess_kurtosis_unbiased <- function(x) {
  x <- as.vector(x)
  T <- length(x)
  x <- x - mean(x)
  #excess_kurt <- mean(x^4)/(mean(x^2))^2 - 3
  excess_kurt <- (T-1)*(T+1)/((T-2)*(T-3)) * (mean(x^4)/(mean(x^2))^2 - 3*(T-1)/(T+1))  # this is better
  excess_kurt_unbiased <- (T-1) / (T-2) / (T-3) * ((T+1)*excess_kurt + 6)  # is this bias correction still necessary?
  return(excess_kurt_unbiased)
}


nu_from_kurtosis <- function(X) {
  kurt <- apply(X, 2, excess_kurtosis_unbiased)
  kappa <- max(0, mean(kurt)/3)
  nu <- 2 / kappa + 4
  # in case of some funny results
  if (nu < .nu_min) nu <- .nu_min
  if (nu > .nu_max) nu <- .nu_max
  return(nu)
}


# estimate nu via Pareto-tail index
alpha_Pareto_tail_index <- function(X, center = FALSE, method = c("WLS", "WLS-stacked", "MLE", "MLE-unbiased")) {
  if (!is.matrix(X)) stop("X must be a matrix.")

  # center data if necessary
  if (center) {
    mu <- colMeans(X)
    X <- X - matrix(mu, nrow(X), ncol(X), byrow = TRUE)
  }
  T <- nrow(X)
  N <- ncol(X)

  # method
  inv_alpha <- switch(match.arg(method),
                      "MLE"          = mean(apply(abs(X), 2, function(x) mean(log(x/min(x))))),
                      "MLE-unbiased" = mean(apply(abs(X), 2, function(x) T/(T-2)*mean(log(x/min(x))))),
                      "WLS"          = mean(apply(abs(X), 2, function(x) mean(log(x/min(x)))/mean(log(T/(1:T))))),
                      "WLS-stacked" = {
                        absXdivXmin <- apply(abs(X), 2, function(x) x/min(x))
                        mean(log(c(absXdivXmin)))/mean(log((N*T)/(1:(N*T))))
                      },
                      stop("Method to estimate Pareto-tail index unknown."))
  return(1/inv_alpha)
}



# estimate nu via MLE
#' @importFrom stats optimize
nu_mle <- function(X, Xc,
                   method = c("MLE-mv-diagcov-resampled", "MLE-mv-cov", "MLE-mv-diagcov",
                              "MLE-mv-scat", "MLE-mv-diagscat", "MLE-mv-diagscat-resampled",
                              "MLE-uv-var-ave", "MLE-uv-scat-ave", "MLE-uv-var-stacked", "test"),
                   Sigma_cov = NULL, Sigma_scatter = NULL) {
  # center data if necessary
  if (missing(Xc)) {
    if (missing(X)) stop("Either X or Xc must be passed.")
    if (!is.matrix(X)) stop("X must be a matrix.")
    mu <- colMeans(X)
    Xc <- X - matrix(mu, nrow(X), ncol(X), byrow = TRUE)
  } else if (!is.matrix(Xc)) stop("Xc must be a matrix.")
  T <- nrow(Xc)
  N <- ncol(Xc)

  # method
  nu <- switch(match.arg(method),
               "MLE-mv-cov" = {  # not good with sample mean and SCM
                 if (is.null(Sigma_cov)) Sigma_cov <- cov(X)
                 delta2_cov <- rowSums(Xc * (Xc %*% solve(Sigma_cov)))
                 negLL <- function(nu) (N*T)/2*log((nu-2)/nu) + ((N + nu)/2)*sum(log(nu + nu/(nu-2)*delta2_cov)) -
                   T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-mv-scat" = {
                 if (is.null(Sigma_scatter)) stop("Scatter matrix must be passed.")
                 delta2 <- rowSums(Xc * (Xc %*% solve(Sigma_scatter)))
                 negLL <- function(nu) ((N + nu)/2)*sum(log(nu + delta2)) - T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-mv-diagcov" = {
                 if (!is.null(Sigma_cov)) var <- diag(Sigma_cov)
                 else var <- apply(Xc^2, 2, sum)/(T-1)
                 delta2_var <- Xc^2 / matrix(var, T, N, byrow = TRUE)
                 delta2_cov_diag <- rowSums(delta2_var)  # this amounts to using a diagonal cov matrix
                 negLL <- function(nu) (N*T)/2*log((nu-2)/nu) + ((N + nu)/2)*sum(log(nu + nu/(nu-2)*delta2_cov_diag)) -
                   T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-mv-diagscat" = {
                 if (is.null(Sigma_scatter)) stop("Scatter matrix must be passed.")
                 diagscat <- diag(Sigma_scatter)
                 delta2_diagscat <- Xc^2 / matrix(diagscat, T, N, byrow = TRUE)
                 delta2 <- rowSums(delta2_diagscat)  # this amounts to using a diagonal cov matrix
                 negLL <- function(nu) ((N + nu)/2)*sum(log(nu + delta2)) - T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-mv-diagcov-resampled" = {  # this method is the winner
                 fT_resampling <- 4
                 N_resampling <- round(1.2*N)
                 var <- apply(Xc^2, 2, sum)/(T-1)
                 delta2_var <- Xc^2 / matrix(var, T, N, byrow = TRUE)
                 # delta2_var_resampled <- t(apply(delta2_var[sample(1:T, size = f_resampling*T, replace = TRUE), ], MARGIN = 1,
                 #                           function(x) x[sample(1:N, size = N, replace = TRUE)]))
                 # delta2_var_resampled <- t(apply(delta2_var[rep.int(1:T, times = f_resampling), ], MARGIN = 1,
                 #                           function(x) x[sample(1:N, size = round(1.2*N), replace = TRUE)]))  #replace=FALSE is not as good
                 # delta2_var_resampled <- t(apply(matrix(rep(t(delta2_var), times = fT_resampling), ncol = N, byrow = TRUE), MARGIN = 1,
                 #                           function(x) x[sample(1:N, size = N_resampling, replace = TRUE)]))  #replace=FALSE is not as good
                 # delta2_cov <- rowSums(delta2_var_resampled)
                 delta2_cov <- apply(matrix(rep(t(delta2_var), times = fT_resampling), ncol = N, byrow = TRUE), MARGIN = 1,
                                     function(x) sum(x[sample(1:N, size = N_resampling, replace = TRUE)]))  #replace=FALSE is not as good
                 # delta2_cov <- c(replicate(fT_resampling,
                 #                           apply(delta2_var, MARGIN = 1, function(x) sum(x[sample(1:N, size = N_resampling, replace = TRUE)]))))
                 T <- T*fT_resampling
                 #N <- N_resampling
                 negLL <- function(nu) (N*T)/2*log((nu-2)/nu) + ((N + nu)/2)*sum(log(nu + nu/(nu-2)*delta2_cov)) -
                   T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-mv-diagscat-resampled" = {
                 if (is.null(Sigma_scatter)) stop("Scatter matrix must be passed.")
                 f_resampling <- 2
                 diagscat <- diag(Sigma_scatter)
                 delta2_diagscat <- Xc^2 / matrix(diagscat, T, N, byrow = TRUE)  # this amounts to using a diagonal cov matrix
                 delta2_diagscat_resampled <- t(apply(delta2_diagscat[rep(1:T, times = f_resampling), ], 1,
                                                      function(x) x[sample(1:N, size = round(1.2*N), replace = TRUE)]))
                 delta2_diagscat_resampled <- delta2_diagscat
                 T <- nrow(delta2_diagscat_resampled)
                 N <- ncol(delta2_diagscat_resampled)
                 delta2 <- rowSums(delta2_diagscat_resampled)
                 negLL <- function(nu) ((N + nu)/2)*sum(log(nu + delta2)) - T*lgamma((N + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               "MLE-uv-var-ave" = {  # not so good with sample mean and SCM
                 var <- apply(Xc^2, 2, sum)/(T-1)
                 delta2_var <- Xc^2 / matrix(var, T, N, byrow = TRUE)
                 negLL_uv_var_ave <- function(nu, delta2_vari) T/2*log((nu-2)/nu) + ((1 + nu)/2)*sum(log(nu + nu/(nu-2)*delta2_vari)) -
                   T*lgamma((1 + nu)/2) + T*lgamma(nu/2) - (nu/2)*T*log(nu)
                 nu_i <- apply(delta2_var, 2, function(delta2_vari)
                   optimize(negLL_uv_var_ave, interval = c(.nu_min, .nu_max), delta2_vari = delta2_vari)$minimum)
                 mean(nu_i)
               },
               "MLE-uv-scat-ave" = {  # not so good
                 if (is.null(Sigma_scatter)) stop("Scatter matrix must be passed.")
                 delta2 <- Xc^2 / matrix(diag(Sigma_scatter), T, N, byrow = TRUE)
                 negLL_uv_scat_ave <- function(nu, delta2_i) ((1 + nu)/2)*sum(log(nu + delta2_i)) -
                   T*lgamma((1 + nu)/2) + T*lgamma(nu/2) - (nu*T/2)*log(nu)
                 nu_i <- apply(delta2, 2, function(delta2_i)
                   optimize(negLL_uv_scat_ave, interval = c(.nu_min, .nu_max), delta2_i = delta2_i)$minimum)
                 mean(nu_i)
               },
               "MLE-uv-var-stacked" = {  # not so good
                 var <- apply(Xc^2, 2, sum)/(T-1)
                 delta2_var <- Xc^2 / matrix(var, T, N, byrow = TRUE)
                 negLL <- function(nu) (N*T)/2*log((nu-2)/nu) + ((1 + nu)/2)*sum(log(nu + nu/(nu-2)*c(delta2_var))) -
                   N*T*lgamma((1 + nu)/2) + N*T*lgamma(nu/2) - (nu/2)*N*T*log(nu)
                 optimize(negLL, interval = c(.nu_min, .nu_max))$minimum
               },
               stop("Method to estimate nu unknown."))
  return(nu)
}

