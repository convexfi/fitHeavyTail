context("Function \"fit_mvt()\"")
#library(testthat)

# # generate the multivariate Student's t data ----------------
# set.seed(123)
# N <- 5
# T <- 50
# nu <- 6
# X <- mvtnorm::rmvt(n = T, sigma = diag(N), df = nu, delta = rep(0, N))
# colnames(X) <- c(1:N)
# X_xts <- xts::as.xts(X, order.by = as.Date("1975-04-28") + 1:nrow(X))
# save(X, X_xts, file = "X.RData", version = 2, compress = "xz")

load("X.RData")

test_that("error control works", {
  expect_error(fit_mvt(X = median), "\"X\" must be a matrix or coercible to a matrix.")
  expect_error(fit_mvt(X = "HongKong"), "\"X\" only allows numerical or NA values.")
  expect_error(fit_mvt(X = 1),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_mvt(X[1:4, ]),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_mvt(X = X, max_iter = -1), "\"max_iter\" must be greater than 1.")
  expect_error(fit_mvt(X = X, factors = -1), "\"factors\" must be no less than 1 and no more than column number of \"X\".")
  expect_error(fit_mvt(X = X, nu = "lala"), "'arg' should be one of \"kurtosis\", \"MLE-diag\", \"MLE-diag-resampled\", \"iterative\"", fixed = TRUE)
  expect_error(fit_mvt(X = X, nu = 1), "Non-valid value for nu.")
})


test_that("default mode works", {
  # mvt_model_check <- fit_mvt(X)
  # save(mvt_model_check, file = "fitted_mvt_check.RData", version = 2, compress = "xz")
  mvt_model <- fit_mvt(X)

  load("fitted_mvt_check.RData")
  expect_equal(mvt_model[c("mu", "cov", "scatter", "converged")],
               mvt_model_check[c("mu", "cov", "scatter", "converged")], tolerance = 0.1)

  # test for xts
  fitted_xts <- fit_mvt(X_xts)
  expect_identical(mvt_model[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")],
                   fitted_xts[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")])

  # test for vector
  fitted_1colmatrix <- fit_mvt(X[, 1])
  fitted_vector <- fit_mvt(as.vector(X[, 1]))
  expect_identical(fitted_1colmatrix[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")],
                   fitted_vector[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")])
})


test_that("Gaussian case fits", {
  mvt_model <- fit_mvt(X, nu = Inf)
  expect_equal(mvt_model$mu, colMeans(X))
  expect_equal(mvt_model$cov, (nrow(X)-1)/nrow(X) * cov(X))
})


test_that("factor structure constraint on scatter matrix works", {
  # mvt_model_factor_check <- fit_mvt(X, factors = 3)
  # save(mvt_model_factor_check, file = "fitted_mvt_factor_check.RData", version = 2, compress = "xz")
  mvt_model_factor <- fit_mvt(X, factors = 3)
  load("fitted_mvt_factor_check.RData")
  expect_equal(mvt_model_factor[c("mu", "cov", "scatter", "psi", "converged")],
               mvt_model_factor_check[c("mu", "cov", "scatter", "psi", "converged")], tolerance = 0.1)
  expect_equal(mvt_model_factor$B %*% t(mvt_model_factor$B),
               mvt_model_factor_check$B %*% t(mvt_model_factor_check$B), tolerance = 0.2)
})


test_that("X with NAs works", {
  X_wNA <- X
  for (i in 1:5) X_wNA[i, i] <- NA
  # mvt_model_wNA_check <- fit_mvt(X_wNA)
  # save(mvt_model_wNA_check, file = "fitted_mvt_wNA_check.RData", version = 2, compress = "xz")
  mvt_model_wNA <- fit_mvt(X_wNA)
  load("fitted_mvt_wNA_check.RData")
  expect_equal(mvt_model_wNA[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")],
               mvt_model_wNA_check[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")])
})


test_that("fixed nu works", {
  # mvt_model_fixednu_check <- fit_mvt(X, nu = fitHeavyTail:::nu_from_kurtosis(X))
  # save(mvt_model_fixednu_check, file = "fitted_mvt_fixednu_kurtosis_check.RData", version = 2, compress = "xz")
  load("fitted_mvt_fixednu_kurtosis_check.RData")
  mvt_model_fixednu <- fit_mvt(X, nu = "kurtosis")
  expect_equal(mvt_model_fixednu[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")],
               mvt_model_fixednu_check[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")])

  # mvt_model_fixednu_check <- fit_mvt(X, nu = fitHeavyTail:::nu_mle(X))
  # save(mvt_model_fixednu_check, file = "fitted_mvt_fixednu_mle_check.RData", version = 2, compress = "xz")
  load("fitted_mvt_fixednu_mle_check.RData")
  mvt_model_fixednu <- fit_mvt(X, nu = "MLE-diag-resampled")
  expect_equal(mvt_model_fixednu[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")],
               mvt_model_fixednu_check[c("mu", "cov", "scatter", "nu", "converged", "num_iterations")], tolerance = 0.3)
})

