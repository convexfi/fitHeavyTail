library(covHeavyTail)
library(tictoc)

set.seed(234)

N <- 20
T <- 4*N
factors <- 5
beta <- matrix(rnorm(N*factors, 10, 1), N, factors)
psi <- rexp(N, 1/10)
Sigma_scale <- beta %*% t(beta) + diag(psi)  # scale matrix (not covariance)
mu <- rnorm(N, 0, 1)
nu <- 7

X <- mvtnorm::rmvt(n = T, sigma = Sigma_scale, df = nu, delta = mu)

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



tic()
fit_old <- momentsStudentt(X)
toc()

# Naming:
# fit_mvt <- momentsStudentt
# fit_mvtFA <- covTFA
# fit_mvskewt


tic()
fit_nom <- fit_mvt(X, ftol = 1e6, return_iterates = TRUE)  # use return_convergence? return_iterations? Or debug?
toc()


# fit_wFA <- covTFA(X, factors = factors, ptol = ptol, ftol = ftol, return_iterates = TRUE)  # check better name for factors?
# fit_wFA_wNA <- covTFA(X_wNA, factors = factors, ptol = ptol, ftol = ftol, return_iterates = TRUE)

fit_old$nu
fit_nom$nu

norm(fit_old$mu - mu, "2") / norm(mu, "2")
norm(fit_nom$mu - mu, "2") / norm(mu, "2")
norm(fit_old$mu - fit_nom$mu, "2") / norm(fit_nom$mu, "2")
norm(fit_old$cov * (fit_old$nu - 2) / fit_old$nu - fit_nom$Sigma, "F") / norm(fit_nom$Sigma, "F")


plot(sapply(fit_nom$iterations_record, function(x) x$nu))
sapply(fit_wFA$iterations_record, function(x) x$nu)
sapply(fit_wFA_wNA$iterations_record, function(x) x$nu)
