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
#'          it can be estimated by maximizing the log-likelihood or a surrogate function (methods \code{"ECME"} and
#'          \code{"ECM"}, respectively), and it can be estimated by maximizing the log-likelihood regularized with a
#'          target \code{nu_target} with weight \code{nu_regcoef > 0} (the regularization term is
#'          \code{nu_regcoef * (nu - nu_target)^2}).
#'
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{cov}: default is the data sample covariance matrix,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         \item{\code{nu}: default is \code{4},}
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
#' @param nu Degrees of freedom (\code{>2}) of the \eqn{t} distribution. If a number is passed,
#'           then \code{nu} will be fixed to this number and will not be further optimized;
#'           if \code{"kurtosis"} is passed, then \code{nu} will be computed from the marginal
#'           kurtosis of the time series; otherwise (default option), \code{nu} will be
#'           estimated via the EM method.
#' @param method_nu String indicating the method for estimating \code{nu} (in case \code{nu} was not passed):
#'                  \itemize{\item{\code{"ECM"}: maximize the Q function w.r.t. \code{nu}}
#'                           \item{\code{"ECME"}: maximize the L function w.r.t. \code{nu}.}}
#'                  This argument is used only when there are no NAs in the data and no factor model is chosen.
#' @param nu_target Number (\code{>=2}) indicating the target for the regularization term for \code{nu}
#'                  in case it is estimated (by default it is obtained via the marginal kurtosis).
#' @param nu_regcoef Number (\code{>=0}) indicating the weight of the regularization term for \code{nu}
#'                   in case it is estimated (default is \code{0}, so no regularion is used).
#' @param return_iterates Logical value indicating whether to record the values of the parameters \code{mu},
#'                        \code{scatter}, and \code{nu} (and possibly the log-likelihood if \code{ftol} is used)
#'                        at each iteration (default is \code{FALSE}).
#' @param verbose Logical value indicating whether to allow the function to print messages (default is \code{FALSE}).
#'
#' @return A list containing possibly the following elements:
#'         \item{\code{mu}}{Mean vector estimate.}
#'         \item{\code{cov}}{Covariance matrix estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{B}}{Factor model loading matrix estimate according to \code{cov = (B \%*\% t(B) + diag(psi)}
#'                         (only if factor model requested).}
#'         \item{\code{psi}}{Factor model idiosynchratic variances estimates according to \code{cov = (B \%*\% t(B) + diag(psi)}
#'                           (only if factor model requested).}
#'         \item{\code{log_likelihood}}{Value of log-likelihood after converge of the estimation algorithm
#'                                      (only if \code{ftol < Inf}).}
#'         \item{\code{iterates_record}}{Iterates of the parameters (\code{mu}, \code{scatter}, \code{nu},
#'                                       and possibly \code{log_likelihood} (if \code{ftol < Inf})) along the iterations
#'                                       (only if \code{return_iterates = TRUE}).}
#'         \item{\code{converged}}{Boolean denoting whether the algorithm has converged (\code{TRUE}) or the maximum number
#'                                 of iterations \code{max_iter} has reached (\code{FALSE}).}
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
fit_mvt <- function(X, initial = NULL, factors = ncol(X),
                    max_iter = 100, ptol = 1e-3, ftol = Inf,
                    nu = NULL, method_nu = c("ECM", "ECME"), nu_target = NULL, nu_regcoef = 0,
                    return_iterates = FALSE, verbose = FALSE) {
  ####### error control ########
  X <- try(as.matrix(X), silent = TRUE)
  if (!is.matrix(X)) stop("\"X\" must be a matrix or coercible to a matrix.")
  if (is.null(colnames(X))) colnames(X) <- paste0("Var", 1:ncol(X))
  if (!all(is.na(X) | is.numeric(X))) stop("\"X\" only allows numerical or NA values.")
  if (ncol(X) <= 1) X <- X[!is.na(X), , drop = FALSE]
  if (nrow(X) == 1) stop("Only T=1 sample!!")
  if (nrow(X) < ncol(X)) stop("Cannot deal with T < N, too few samples.")
  factors <- round(factors)
  max_iter <- round(max_iter)
  if (factors < 1 || factors > ncol(X)) stop("\"factors\" must be no less than 1 and no more than column number of \"X\".")
  if (max_iter < 1) stop("\"max_iter\" must be greater than 1.")
  ##############################

  T <- nrow(X)
  N <- ncol(X)
  method_nu <- match.arg(method_nu)
  X_has_NA <- anyNA(X)
  FA_struct <- (factors != N)
  optimize_nu <- is.null(nu)
  if (!optimize_nu) {
    if (nu == Inf) nu <- 1e15  # for numerical stability (for the Gaussian case)
    else if (nu == "kurtosis") {  # estimate nu if argument nu = "kurtosis"
      nu <- est_nu_kurtosis(X)
      if (verbose) message(sprintf("Automatically set nu = %.2f", nu))
    } else if (!is.numeric(nu) || nu <=2 ) stop("Non-valid value for nu.")
  } else  # choose nu_target if necessary
    if (is.null(nu_target)) {
      if (nu_regcoef > 0) {  # really need nu_target
        nu_target <- est_nu_kurtosis(X)
        if (verbose) message(sprintf("Automatically choose a target nu = %.2f", nu_target))
      } else  # no need to have nu_target, but assign a value to simplify following codes
        nu_target <- 0
    } else if (nu_target < 2) stop("Non-valid value for nu_target.")

  # initialize all parameters
  alpha <- 1  # an extra variable for PX-EM acceleration
  if (optimize_nu) nu <- if (is.null(initial$nu)) 4 else initial$nu
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

  # aux function to save iterates
  snapshot <- function() {
    if (ftol < Inf) list(mu = mu, scatter = Sigma, nu = nu, log_likelihood = log_likelihood)
    else list(mu = mu, scatter = Sigma, nu = nu)
  }


  # loop
  if (return_iterates) iterates_record <- list(snapshot())
  for (iter in 1:max_iter) {
    # record the current status
    Sigma_old <- Sigma
    mu_old <- mu
    nu_old <- nu
    if (ftol < Inf) log_likelihood_old <- log_likelihood

    ## -------------- E-step --------------
    if (X_has_NA || FA_struct)
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
      Q_nu <- function(nu) { - (nu/2)*log(nu/2) + lgamma(nu/2) - (nu/2)*(Q$ave_E_logtau - Q$ave_E_tau) + nu_regcoef * (nu/(nu-2) - nu_target/(nu_target-2))^2 }
      if (optimize_nu) nu <- optimize(Q_nu, interval = c(2 + 1e-12, 100))$minimum
    } else {
      mu <- ave_E_tau_X / ave_E_tau
      alpha <- ave_E_tau  # acceleration
      X_ <- X - matrix(mu, T, N, byrow = TRUE)  # this is slower: sweep(X, 2, FUN = "-", STATS = mu)  #X_ <- X - rep(mu, each = TRUE)  # this is wrong?
      ave_E_tau_XX <- (1/T) * crossprod(sqrt(E_tau) * X_)  # (1/T) * t(X_) %*% diag(E_tau) %*% X_
      Sigma <- ave_E_tau_XX / alpha
      if (optimize_nu)
        nu <- switch(method_nu,
                     "ECM" = {  # based on minus the Q function of nu
                       E_log_tau_minus_E_tau <- T*(digamma((N+nu)/2) - log((N+nu)/2)) + sum(log(E_tau) - E_tau)
                       Q_nu <- function(nu) { - T*(nu/2)*log(nu/2) + T*lgamma(nu/2) - (nu/2)*E_log_tau_minus_E_tau +
                                              nu_regcoef * (nu/(nu-2) - nu_target/(nu_target-2))^2 }
                       optimize(Q_nu, interval = c(2 + 1e-12, 100))$minimum
                       },
                     "ECME" = {  # based on minus log-likelihood of nu with mu and sigma fixed to mu[k+1] and sigma[k+1]
                       tmp <- rowSums(X_ * (X_ %*% inv(Sigma)))  # diag( X_ %*% inv(Sigma) %*% t(X_) )
                       LL_nu <- function(nu) { - sum(- ((nu+N)/2)*log(nu+tmp) + lgamma( (nu+N)/2 ) - lgamma(nu/2) + (nu/2)*log(nu)) +
                                               nu_regcoef * (nu/(nu-2) - nu_target/(nu_target-2))^2 }
                       optimize(LL_nu, interval = c(2 + 1e-12, 100))$minimum
                       },
                     stop("Method unknown."))
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
  if (verbose) message(sprintf("Number of iterations for mvt estimator = %d\n", iter))

  ## -------- return variables --------
  #Sigma <- T/(T-1) * Sigma  # unbiased estimator
  vars_to_be_returned <- list("mu"          = mu,
                              "cov"         = nu/(nu-2) * Sigma,
                              "scatter"     = Sigma,
                              "nu"          = nu)
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
  vars_to_be_returned$converged <- (iter < max_iter)
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
  excess_kurt <- (T+1)*(T-1) * ((sum(x^4)/T)/(sum(x^2)/T)^2 - 3*(T-1)/(T+1))/((T-2)*(T-3))
  excess_kurt_unbiased <- (T-1) / (T-2) / (T-3) * ((T+1)*excess_kurt + 6)  # is this bias correction still necessary?
  return(excess_kurt_unbiased)
}


est_nu_kurtosis <- function(X) {
  kurt <- apply(X, 2, excess_kurtosis_unbiased)
  kappa <- max(0, mean(kurt)/3)
  nu <- 2 / kappa + 4
  # in case of some funny results
  if (nu < 2) nu <- 2
  if (nu > 100) nu <- 100
  return(nu)
}
