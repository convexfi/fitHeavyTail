context("Function \"fit_mvst()\"")
#library(testthat)

# # generate the multivariate skew t data ----------------
# library(mvtnorm)
# set.seed(123)
# N <- 5
# T <- 1000
# nu <- 6
# mu <- rnorm(N)
# scatter <- diag(N)
# gamma <- rnorm(N)   # skewness vector
# taus <- rgamma(n = T, shape = nu/2, rate = nu/2)
# X <- matrix(data = mu, nrow = T, ncol = N, byrow = TRUE) +
#      matrix(data = gamma, nrow = T, ncol = N, byrow = TRUE) / taus +
#      rmvnorm(n = T, mean = rep(0, N), sigma = scatter) / sqrt(taus)
# colnames(X) <- c(1:N)
# X_xts <- xts::as.xts(X, order.by = as.Date("1975-04-28") + 1:nrow(X))
# save(X, X_xts, nu, mu, scatter, gamma, file = "X_mvst.RData", version = 2, compress = "xz")


load("X_mvst.RData")


test_that("error control works", {
  expect_error(fit_mvt(X = median), "\"X\" must be a matrix or coercible to a matrix.")
  expect_error(fit_mvt(X = "HongKong"), "\"X\" only allows numerical or NA values.")
  expect_error(fit_mvt(X = 1),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_mvt(X[1:4, ]),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_mvt(X = X, max_iter = -1), "\"max_iter\" must be greater than 1.")
})


test_that("default mode works", {
  # mvst_model_check <- fit_mvst(X, ptol = 1e-4)
  # save(mvst_model_check, file = "fitted_mvst_check.RData", version = 2, compress = "xz")
  load("fitted_mvst_check.RData")

  mvst_model <- fit_mvst(X, ptol = 1e-4)

  expect_equal(mvst_model[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               mvst_model_check[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               tolerance = 1e-5)

  # test for xts
  mvst_model <- fit_mvst(X)
  fitted_xts <- fit_mvst(X_xts)
  expect_equal(mvst_model[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               fitted_xts[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               tolerance = 1e-5)

  # test for vector
  fitted_1colmatrix <- fit_mvst(X[, 1])
  fitted_vector <- fit_mvst(as.vector(X[, 1]))
  expect_equal(fitted_1colmatrix[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               fitted_vector[c("mu", "gamma", "scatter", "nu", "mean", "cov", "converged", "num_iterations")],
               tolerance = 1e-5)


  # mvst_model <- fit_mvst(X, ftol = 1e-6, return_iterates = TRUE)
  # plot(sapply(mvst_model$iterates_record, function(x) x$obj))
})


test_that("Gaussian case fits", {
  # mvt_model <- fit_mvt(X, nu = Inf)
  # expect_equal(mvt_model$mu, colMeans(X))
  # expect_equal(mvt_model$cov, (nrow(X)-1)/nrow(X) * cov(X))

  mvt_model  <- fit_mvt(X,  nu = 1e5, ftol = 1e-3)
  mvst_model <- fit_mvst(X, nu = 1e5, gamma = rep(0, ncol(X)), ftol = 1e-3)
  expect_equal(mvt_model$mu,  mvst_model$mu)
  expect_equal(mvt_model$cov, mvst_model$cov)


  # sanity check for large nu:
  mvt_model  <- fit_mvt(X,  nu = 10, ftol = 1e-3)
  mvst_model <- fit_mvst(X, nu = 10, gamma = rep(0, ncol(X)), ftol = 1e-3)
  expect_equal(mvt_model$mu,  mvst_model$mu)
  expect_equal(mvt_model$cov, mvst_model$cov)


  expect_silent(
    mvst_model <- fit_mvst(X, nu = 1e5, gamma = rep(1e-6, ncol(X)), ftol = 1e-3)
  )
})



test_that("bounds on nu work", {
  options(nu_min = 2.5)
  mvst_model <- fit_mvst(X)
  expect_false(mvst_model$nu > 9.37)

  options(nu_min = 9.37)
  mvst_model <- fit_mvst(X)
  expect_true(mvst_model$nu > 9.37 - 1e-8)

  options(nu_min = 2.5)
})



test_that("Bessel functions work", {
  res1 <- sum(dmvst_orig(X = X, nu = 5.7, gamma = gamma, mu = mu, scatter = scatter))
  res2 <- sum(dmvst(X = X, nu = 5.7, gamma = gamma, mu = mu, scatter = scatter))
  expect_equal(res1, res2)

  res1 <- sum(dmvst_orig(X = X, nu = 15.7, gamma = gamma, mu = mu, scatter = scatter))
  res2 <- sum(dmvst(X = X, nu = 15.7, gamma = gamma, mu = mu, scatter = scatter))
  expect_equal(res1, res2)

  res1 <- sum(dmvst_orig(X = X, nu = 60.7, gamma = gamma, mu = mu, scatter = scatter))
  res2 <- sum(dmvst(X = X, nu = 60.7, gamma = gamma, mu = mu, scatter = scatter))
  expect_equal(res1, res2)
  res1 <- besselK_ratio(x = 10, nu = 70)
  res2 <- besselK(x=10, nu = 71)/besselK(x=10, nu = 70)
  expect_equal(res1, res2)
  res1 <- besselK_ratio(x = 0.01, nu = 5)
  res2 <- besselK(x=0.01, nu = 6)/besselK(x=0.01, nu = 5)
  expect_equal(res1, res2)
  res1 <- besselK_ratio(x = 10, nu = 5)
  res2 <- besselK(x=10, nu = 6)/besselK(x=10, nu = 5)
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 0.01, nu = 5)
  res2 <- log(besselK(0.01, nu = 5))
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 0.01, nu = -5)
  res2 <- log(besselK(0.01, nu = -5))
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 10, nu = 5)
  res2 <- log(besselK(10, nu = 5))
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 10, nu = -5)
  res2 <- log(besselK(10, nu = -5))
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 10, nu = 70)
  res2 <- log(besselK(10, nu = 70))
  expect_equal(res1, res2)
  res1 <- log_besselK(x = 10, nu = -70)
  res2 <- log(besselK(10, nu = -70))
  expect_equal(res1, res2)
})
