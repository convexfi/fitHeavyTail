context("Function \"fit_Tyler()\"")
#library(testthat)


load("X_mvt.RData")
# recall:
# scatter_true <- diag(5)
# nu_true <- 6
# cov_true <- nu_true/(nu_true-2) * scatter_true


test_that("error control works", {
  expect_error(fit_Tyler(X = median),
               "\"X\" must be a matrix or coercible to a matrix.")
  expect_error(fit_Tyler(X = "HongKong"),
               "\"X\" only allows numerical or NA values.")
  expect_error(fit_Tyler(X = 1),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_Tyler(X[1:4, ]),
               "Cannot deal with T <= N (after removing NAs), too few samples.", fixed = TRUE)
  expect_error(fit_Tyler(X = X, max_iter = -1),
               "\"max_iter\" must be greater than 1.")
})


test_that("cov estimate works", {
  # test against fit_mvt()
  fitted_Tyler <- fit_Tyler(X)
  fitted_mvt   <- fit_mvt(X)
  expect_equal(fitted_Tyler$cov, fitted_mvt$cov, tolerance = 0.3)
  # norm(fitted_Tyler$cov - cov_true, "F")
  # norm(fitted_mvt$cov - cov_true, "F")

  # # plotting convergence
  # fitted_Tyler <- fit_Tyler(X, ftol = 1, verbose = TRUE, return_iterates = TRUE)
  # fitHeavyTail:::plotConvergence(fitted_Tyler)


  # test against saved results
  # fitted_Tyler_check <- fit_Tyler(X)
  # save(fitted_Tyler_check, file = "fitted_Tyler_check.RData", version = 2, compress = "xz")
  load("fitted_Tyler_check.RData")
  expect_equal(fitted_Tyler[c("mu", "cov", "converged")], fitted_Tyler_check[c("mu", "cov", "converged")])

  # test for xts
  fitted_xts <- fit_Tyler(X_xts)
  expect_identical(fitted_Tyler[c("mu", "cov", "converged", "num_iterations")],
                   fitted_xts[c("mu", "cov", "converged", "num_iterations")])

  # test for vector
  fitted_1colmatrix <- fit_Tyler(X[, 1])
  fitted_vector <- fit_Tyler(as.vector(X[, 1]))
  expect_identical(fitted_1colmatrix[c("mu", "cov", "converged", "num_iterations")],
                   fitted_vector[c("mu", "cov", "converged", "num_iterations")])
})


test_that("X with NAs works", {
  X_wNA <- X
  for (i in 1:5) X_wNA[i, i] <- NA

  fitted_Tyler <- fit_Tyler(X[-c(1:5), ])
  expect_message(fitted_Tyler_wNA <- fit_Tyler(X_wNA, verbose = TRUE), "X contains NAs, dropping those observations.")
  expect_identical(fitted_Tyler[c("mu", "cov", "converged")], fitted_Tyler_wNA[c("mu", "cov", "converged")])
})

