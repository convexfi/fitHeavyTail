library(covHeavyTail)

set.seed(234)

p <- 100
n <- p * 4
r <- 5
beta <- matrix(rnorm(p*r, 10, 1), p, r)
psi <- rexp(p, 1/10)
sigma <- beta %*% t(beta) + diag(psi, p)
mu <- rnorm(p, 0, 1)
nu <- 7

X <- mvtnorm::rmvt(n = n, sigma = sigma, df = nu, delta = mu)

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

max_iter <- 100
ptol <- 0
ftol <- 1e-6

fit_old <- momentsStudentt(X)

fit_nom <- covTFA(X, ptol = ptol, ftol = ftol, procedure = TRUE)
fit_wFA <- covTFA(X, r = r, ptol = ptol, ftol = ftol, procedure = TRUE)
fit_wFA_wNA <- covTFA(X_wNA, r = r, ptol = ptol, ftol = ftol, procedure = TRUE)

fit_old$nu
fit_nom$nu

norm(fit_old$mu - fit_nom$mu, "2") / norm(fit_nom$mu, "2")
norm(fit_old$cov * (fit_old$nu - 2) / fit_old$nu - fit_nom$Sigma, "F") / norm(fit_nom$Sigma, "F")


sapply(fit_nom$proc, function(x) x$nu)
sapply(fit_wFA$proc, function(x) x$nu)
sapply(fit_wFA_wNA$proc, function(x) x$nu)
