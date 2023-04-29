
# estimate nu via the kurtosis of each variable
excess_kurtosis_unbiased <- function(x) {
  x <- as.vector(x)
  T <- length(x)
  x <- x - mean(x)
  #excess_kurt <- mean(x^4)/(mean(x^2))^2 - 3   # vanilla estimator
  excess_kurt <- (T-1)*(T+1)/((T-2)*(T-3)) * (mean(x^4)/(mean(x^2))^2 - 3*(T-1)/(T+1))  # including corrections from RMT
  excess_kurt_unbiased <- (T-1) / (T-2) / (T-3) * ((T+1)*excess_kurt + 6)  # bias correction again from RMT
  return(excess_kurt_unbiased)
}


nu_from_average_marginal_kurtosis <- function(X, remove_outliers = TRUE) {
  X <- na.omit(as.matrix(X))
  N <- ncol(X)

  kurt_columns  <- apply(X, 2, excess_kurtosis_unbiased)
  kappa_columns <- pmax(-2/(N+2)*0.99, kurt_columns/3)

  # eliminate outliers
  if (remove_outliers) {
    nu_columns <- 2/max(1e-10, kappa_columns) + 4
    idx <- nu_columns >= getOption("nu_min") & nu_columns <= getOption("nu_max")
    if (length(idx) > 0)
      kappa_columns <- kappa_columns[idx]
  }

  ave_kappa <- max(1e-10, mean(kappa_columns, na.rm = TRUE))
  nu <- 2/ave_kappa + 4
  min(getOption("nu_max"), max(getOption("nu_min"), nu))
}


# nu_from_average_marginal_kurtosis_whitened <- function(X) {
#   X <- na.omit(as.matrix(X))
#   N <- ncol(X)
#   #Sigma12 <- chol(cov(X))    # norm(cov(X) - t(Sigma12) %*% Sigma12, "F")
#   #X_ <- t(Sigma12 %*% t(X))
#   #X_ <- t(solve(t(Sigma12), t(X)))
#   #X_ <- t(forwardsolve(t(Sigma12), t(X)))
#
#   # eig <- eigen((1 - 0.1)*cov(X) + 0.1*diag(N))
#   # X_ <- t(t(eig$vectors) %*% t(X))
#
#   X <- scale(X)
#   Sigma12 <- chol(cov(X))    # norm(cov(X) - t(Sigma12) %*% Sigma12, "F")
#   X_ <- t(Sigma12 %*% t(X))
#
#   nu_from_average_marginal_kurtosis(X_)
# }


#
# This one is not better than the vanilla kurtosis method
#
# nu_from_average_marginal_kurtosis_resampled <- function(X) {
#   X <- na.omit(as.matrix(X))
#   T <- nrow(X)
#   N <- ncol(X)
#   x <- as.vector(scale(X))
#
#   X_ <- matrix(NA, 2*T, 2*N)
#   for (i in 1:ncol(X_))
#     X_[, i] <- x[sample(length(x), nrow(X_), replace = FALSE)]
#
#   nu_from_average_marginal_kurtosis(X_)
# }



#
# This one is not better than: nu_from_average_marginal_kurtosis_whitened()
#
# nu_from_average_marginal_kurtosis_virtual_vars <- function(X) {
#   X <- na.omit(as.matrix(X))
#   N <- ncol(X)
#   X <- scale(X)
#   num_virtual_observations <- 50*N
#
#   # generate random directions
#   Omega <- t(rmvnorm(num_virtual_observations, sigma = diag(N)))
#
#   # augment observations in X
#   X_augmented <- cbind(X, X %*% Omega)
#
#   nu_from_average_marginal_kurtosis(X_augmented)
# }


nu_from_cross_cumulants <- function(X) {
  X <- na.omit(as.matrix(X))
  N <- ncol(X)
  T <- nrow(X)

  X <- scale(X)
  X2 <- X^2
  sigma2 <- colMeans(X2)
  S <- (1/T)*crossprod(X)
  S2 <- (1/T)*crossprod(X2)
  kappas_offdiag <- c()
  for (i in 1:(N-1))
    kappas_offdiag <- c(kappas_offdiag,
                        S2[i, (i+1):N] / (sigma2[i] * sigma2[(i+1):N] + 2*S[i, (i+1):N]^2) - 1)

  kappas_offdiag <- pmax(-2/(N+2)*0.99, kappas_offdiag)
  ave_kappa <- max(1e-10, mean(kappas_offdiag))
  nu <- 2/ave_kappa + 4
  min(getOption("nu_max"), max(getOption("nu_min"), nu))
}


nu_from_all_cumulants <- function(X) {
  nu_marginals <- nu_from_average_marginal_kurtosis(X)
  nu_cumulants <- nu_from_cross_cumulants(X)

  (nu_marginals + nu_cumulants)/2
}




# [AshurbekovaCarleveForbesAchard2020]
nu_Hill_estimator <- function(X) {
  T <- nrow(X)
  nu_max <- 12
  b <- 4/(nu_max + 4)
  kn <- floor(T^b)
  norm_X <- sqrt(rowSums(X^2))
  sorted_norm_X <- sort(norm_X, decreasing = TRUE)
  nu <- 1/mean(log(sorted_norm_X[1:kn]/sorted_norm_X[kn+1]))
  min(getOption("nu_max"), max(getOption("nu_min"), nu))
}


# estimate nu via Pareto-tail index (this is related to Hill estimator)
nu_Pareto_tail_index <- function(X, center = FALSE, method = c("WLS", "WLS-stacked", "MLE", "MLE-unbiased")) {
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

  min(getOption("nu_max"), max(getOption("nu_min"), 1/inv_alpha))
}
