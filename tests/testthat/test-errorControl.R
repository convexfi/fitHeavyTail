context("Checking package error control")

load("X.RData")

test_that("Error control test for function \"fit_mvt()\"", {

  expect_error(fit_mvt(X = median), "\"X\" must be a matrix or can be converted to a matrix.")

  expect_error(fit_mvt(X = "HongKong"), "\"X\" only allows numerical or NA values.")

  expect_error(fit_mvt(X = 1), "Only T=1 sample!!")

  expect_error(fit_mvt(X = X, factors = -1), "\"factors\" must be no less than 1 and no more than column number of \"X\".")

  expect_error(fit_mvt(X = X, max_iter = -1), "\"max_iter\" must be greater than 1.")

})
