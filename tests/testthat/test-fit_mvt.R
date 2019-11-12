context("Checking function \"fit_mvt()\"")

# generate the multivariate Student's t data ----------------
# set.seed(123)
# N <- 10
# T <- 50
# nu <- 6
# X <- mvtnorm::rmvt(n = T, sigma = diag(N), df = nu, delta = rep(0, N))
# save(X, file = "tests/testthat/X.RData", version = 2, compress = "xz")

load("X.RData")

# default mode --------------------------------------------------------------------------------------------------------
test_that("test the default mode", {
  # mvt_model_check <- fit_mvt(X)
  # save(mvt_model_check, file = "tests/testthat/mvt_model_check.RData", version = 2, compress = "xz")
  mvt_model <- fit_mvt(X)
  load("mvt_model_check.RData")
  expect_equal(mvt_model, mvt_model_check)
})

# factor structure constraint -----------------------------------------------------------------------------------------
test_that("with factor structure constraint on scatter matrix", {
  # mvt_model_factor_check <- fit_mvt(X, 5)
  # save(mvt_model_factor_check, file = "tests/testthat/mvt_model_factor_check.RData", version = 2, compress = "xz")
  mvt_model_factor <- fit_mvt(X, 5)
  load("mvt_model_factor_check.RData")
  expect_equal(mvt_model_factor, mvt_model_factor_check)
})

# missing data in X ---------------------------------------------------------------------------------------------------
test_that("X contains missing values", {
  X_wNA <- X
  for (i in 1:5) X_wNA[i, i] <- NA
  # mvt_model_missing_check <- fit_mvt(X_wNA)
  # save(mvt_model_missing_check, file = "tests/testthat/mvt_model_missing_check.RData", version = 2, compress = "xz")
  mvt_model_missing <- fit_mvt(X_wNA)
  load("mvt_model_missing_check.RData")
  expect_equal(mvt_model_missing, mvt_model_missing_check)
})

# fix nu --------------------------------------------------------------------------------------------------------------
test_that("fix nu", {
  # nu_kurtosis <- fitHeavyTail:::est_nu_kurtosis(X)
  # mvt_model_fixnu_check <- fit_mvt(X, nu = nu_kurtosis)
  # save(mvt_model_fixnu_check, file = "tests/testthat/mvt_model_fixnu_check.RData", version = 2, compress = "xz")
  load("mvt_model_fixnu_check.RData")
  mvt_model_fixnu <- fit_mvt(X, nu = "kurtosis")
  expect_equal(mvt_model_fixnu, mvt_model_fixnu_check)
})

# regularize nu -------------------------------------------------------------------------------------------------------
test_that("regularize nu", {
  # nu_kurtosis <- fitHeavyTail:::est_nu_kurtosis(X)
  # mvt_model_regnu_check <- fit_mvt(X, nu_target = nu_kurtosis, nu_regcoef = 1)
  # save(mvt_model_regnu_check, file = "tests/testthat/mvt_model_regnu_check.RData", version = 2, compress = "xz")
  load("mvt_model_regnu_check.RData")
  mvt_model_regnu <- fit_mvt(X, nu_regcoef = 1)
  expect_equal(mvt_model_regnu, mvt_model_regnu_check)
})
