library(fitHeavyTail)
library(tictoc)

set.seed(234)

N <- 100
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

method <- "ECM"
fit_old <- fit_mvt(X, ptol = 0, ftol = 0, return_iterates = TRUE)
fit_nom <- fit_mvt(X, ptol = 0, ftol = 0, return_iterates = TRUE)

fit_wFA <- fit_mvt(X, factors = factors, ftol = 1e6, return_iterates = TRUE)
fit_wFA_wNA <- fit_mvt(X_wNA, factors = factors, ftol = 1e6, return_iterates = TRUE)
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


# estimate nu from subspace of X ---------------------------
library(mvtnorm)
library(fitHeavyTail)
library(ggplot2)

N <- 50
T <- 1.1*N
nu <- 5
mu <- rep(0, N)

# generate data
set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")

X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data

res_mvt_noreg <- fit_mvt(X, return_iterates = TRUE)
res_mvt_wreg  <- fit_mvt(X, nu_regcoef = 2e-2, return_iterates = TRUE)

p_mvt_noreg <- fitHeavyTail:::plotConvergence(res_mvt_noreg)
do.call(gridExtra::grid.arrange, c(p_mvt_noreg, ncol = 1))

p_mvt_wreg <- fitHeavyTail:::plotConvergence(res_mvt_wreg)
do.call(gridExtra::grid.arrange, c(p_mvt_wreg, ncol = 1))
