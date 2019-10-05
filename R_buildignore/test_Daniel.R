library(mvtnorm)
library(fitHeavyTail)
library(ggplot2)

N <- 20
T <- 1.1*N
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


res_Studentt_nu6 <- fit_mvt(X, nu = 6, method = "ECME", ftol = 1e6, return_iterates = TRUE)
res_Studentt <- fit_mvt(X, method = "ECME", ftol = 1e6, return_iterates = TRUE)

p_Studentt_nu6 <- fitHeavyTail:::plotConvergence(res_Studentt_nu6)
p_Studentt <- fitHeavyTail:::plotConvergence(res_Studentt)

# print everything in one single plot
do.call(gridExtra::grid.arrange, c(p_Studentt_nu6, ncol = 1))
do.call(gridExtra::grid.arrange, c(p_Studentt, ncol = 1))




# #
# # comparison with old function:
# #
# res_old <- momentsStudentt(X, max_iter = 50, verbose = TRUE)
# res_new <- fit_mvt(X, max_iter = 50, ftol = 1e6, return_iterates = TRUE)
#
# # res_old <- momentsStudentt(X, method = "ECME", max_iter = 40, verbose = TRUE)
# # res_new <- fit_mvt(X, method = "ECME", max_iter = 40, ftol = 1e6, return_iterates = TRUE)
#
# # res_old <- momentsStudentt(X, nu = 4, max_iter = 20, verbose = TRUE)
# # res_new <- fit_mvt(X, nu = 4, max_iter = 20, ftol = 1e6, return_iterates = TRUE)
#
# res_old$nu
# res_new$nu
# max(abs(res_old$mu - res_new$mu) / abs(res_old$mu))
# max(abs(res_old$cov - res_new$cov) / abs(res_old$cov))
#
# plot(res_old$obj_value_record, type = "b", col = "blue")
# lines(sapply(res_new$iterations_record, function(x) x$log_likelihood), col = "red")





# #
# # cpu time comparison: the old function is slower because it always computes the log-likelihood and returns its iterations
# #
# library(microbenchmark)
#
# op <- microbenchmark(
#   old = res_old <- momentsStudentt(X, max_iter = 50),
#   new = res_new <- fit_mvt(X, max_iter = 50),
#   times = 100L)
# print(op)
# autoplot(op)








