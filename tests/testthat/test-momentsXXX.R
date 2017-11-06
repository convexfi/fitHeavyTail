data(heavy_data)

test_that("NA check", {
  X_na <- heavy_data$X
  X_na[10, 10] <- NA
  expect_error(momentsStudentt(X_na))
  expect_error(momentsCauchy(X_na))
  expect_error(momentsTyler(X_na))
})


test_that("input dimension check", {
  expect_error(momentsStudentt(heavy_data$X[1, ]))  # T > 1?
  expect_error(momentsStudentt(heavy_data$X[, 1]))  # N > 1?

  expect_error(momentsCauchy(heavy_data$X[1, ]))  # T > 1?
  expect_error(momentsCauchy(heavy_data$X[, 1]))  # N > 1?

  expect_error(momentsTyler(heavy_data$X[1, ]))  # T > 1?
  expect_error(momentsTyler(heavy_data$X[, 1]))  # N > 1?
})


test_that("output dimension check", {
  res <- momentsStudentt(heavy_data$X)
  expect_equal(length(res$mu), heavy_data$N)
  expect_equal(nrow(res$cov), heavy_data$N)
  expect_equal(ncol(res$cov), heavy_data$N)
  expect_equal(length(res$nv), 1)

  res <- momentsCauchy(heavy_data$X)
  expect_equal(length(res$mu), heavy_data$N)
  expect_equal(nrow(res$cov), heavy_data$N)
  expect_equal(ncol(res$cov), heavy_data$N)

  res <- momentsTyler(heavy_data$X)
  expect_equal(length(res$mu), heavy_data$N)
  expect_equal(nrow(res$cov), heavy_data$N)
  expect_equal(ncol(res$cov), heavy_data$N)
})


test_that("mu and cov estimation check", {
  res <- momentsStudentt(heavy_data$X)
  expect_lt(norm(res$mu - heavy_data$mu, "2"), 0.2390831 + 0.001)
  expect_lt(norm(res$cov - heavy_data$cov, "F"), 5.651713 + 0.001)

  res <- momentsCauchy(heavy_data$X)
  expect_lt(norm(res$mu - heavy_data$mu, "2"), 0.2377818 + 0.001)
  expect_lt(norm(res$cov - heavy_data$cov, "F"), 6.574427 + 0.001)

  res <- momentsTyler(heavy_data$X)
  expect_lt(norm(res$mu - heavy_data$mu, "2"), 0.7242409 + 0.001)
  expect_lt(norm(res$cov - heavy_data$cov, "F"), 5.826894 + 0.001)
})


test_that("Gaussian mu and cov estimation check", {
  res <- momentsStudentt(heavy_data$X, nv = 1e15)
  expect_lt(norm(res$mu - colMeans(heavy_data$X), "2"), 1e-12)
  expect_lt(norm(res$cov - ((heavy_data$T-1)/heavy_data$T)*cov(heavy_data$X), "F"), 1e-11)
})
