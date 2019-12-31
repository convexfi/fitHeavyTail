context("Function \"fit_Cauchy()\"")
#library(testthat)


load("X.RData")

test_that("error control works", {
  expect_error(fit_Cauchy(X = median),
               "\"X\" must be a matrix or coercible to a matrix.")
  expect_error(fit_Cauchy(X = "HongKong"),
               "\"X\" only allows numerical or NA values.")
  expect_error(fit_Cauchy(X = 1),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_Cauchy(X[1:4, ]),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_Cauchy(X = X, max_iter = -1),
               "\"max_iter\" must be greater than 1.")
})


test_that("cov estimate works", {
  # test against fit_mvt()
  fitted_Cauchy <- fit_Cauchy(X)
  fitted_mvt <- fit_mvt(X)
  expect_equal(fitted_Cauchy$cov, fitted_mvt$cov, tolerance = 0.4)

  # # plotting convergence
  # fitted_Cauchy <- fit_Cauchy(X, ftol = 1, verbose = TRUE, return_iterates = TRUE)
  # fitHeavyTail:::plotConvergence(fitted_Cauchy)

  # test agains saved results
  # fitted_Cauchy_check <- fit_Cauchy(X)
  # save(fitted_Cauchy_check, file = "fitted_Cauchy_check.RData", version = 2, compress = "xz")
  load("fitted_Cauchy_check.RData")
  # expect_identical(fitted_Cauchy, fitted_Cauchy_check)
  expect_equal(fitted_Cauchy[c("mu", "cov", "scatter","converged")], fitted_Cauchy_check[c("mu", "cov", "scatter","converged")])

  # test for xts
  fitted_xts <- fit_Cauchy(X_xts)
  expect_identical(fitted_Cauchy[c("mu", "cov", "scatter", "converged", "num_iterations")],
                   fitted_xts[c("mu", "cov", "scatter", "converged", "num_iterations")])

  # test for vector
  fitted_1colmatrix <- fit_Cauchy(X[, 1])
  fitted_vector <- fit_Cauchy(as.vector(X[, 1]))
  expect_identical(fitted_1colmatrix[c("mu", "cov", "scatter", "converged", "num_iterations")],
                   fitted_vector[c("mu", "cov", "scatter", "converged", "num_iterations")])
})


test_that("X with NAs works", {
  X_wNA <- X
  for (i in 1:5) X_wNA[i, i] <- NA

  fitted_Cauchy <- fit_Cauchy(X[-c(1:5), ])
  expect_message(fitted_Cauchy_wNA <- fit_Cauchy(X_wNA, verbose = TRUE), "X contains NAs, dropping those observations.")
  expect_identical(fitted_Cauchy[c("mu", "cov", "scatter", "converged")], fitted_Cauchy_wNA[c("mu", "cov", "scatter", "converged")])
})

