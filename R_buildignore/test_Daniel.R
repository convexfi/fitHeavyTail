library(mvtnorm)
library(fitHeavyTail)
library(ggplot2)
library(reshape2)

N <- 5
T <- 1.5*N
nu <- 4
mu <- rep(0, N)

# generate data
set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")

X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
#X <- rmvnorm(n = T, mean = mu, sigma = Sigma)  # Gaussian data


res_Daniel <- momentsStudentt(X, max_iter = 50, verbose = TRUE)
res_Rui <- fit_mvt(X, max_iter = 50, ftol = 1e6, return_iterates = TRUE)

# res_Daniel <- momentsStudentt(X, method = "ECME", max_iter = 40, verbose = TRUE)
# res_Rui <- fit_mvt(X, method = "ECME", max_iter = 40, ftol = 1e6, return_iterates = TRUE)

# res_Daniel <- momentsStudentt(X, nu = 4, max_iter = 20, verbose = TRUE)
# res_Rui <- fit_mvt(X, nu = 4, max_iter = 20, ftol = 1e6, return_iterates = TRUE)


res_Daniel$nu
res_Rui$nu
max(abs(res_Daniel$mu - res_Rui$mu) / abs(res_Daniel$mu))
max(abs(res_Daniel$cov - res_Rui$cov) / abs(res_Daniel$cov))

plot(res_Daniel$obj_value_record, type = "b", col = "blue")
lines(sapply(res_Rui$iterations_record, function(x) x$log_likelihood), col = "red")





#
# cpu time comparison: the old function is slower because it always computes the log-likelihood and returns its iterations
#
library(microbenchmark)

op <- microbenchmark(
  old = res_Daniel <- momentsStudentt(X, max_iter = 50),
  new = res_Rui <- fit_mvt(X, max_iter = 50),
  times = 100L)
print(op)
autoplot(op)








