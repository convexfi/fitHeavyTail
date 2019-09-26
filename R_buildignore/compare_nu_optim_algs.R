# This file is to compare the time consuming when solve the optimal nu
# 1. stats::optimize()
# 2. bisection method
library(tictoc)

tmp <- -1.2

Q_nu <- function(nu) { - (nu/2)*log(nu/2) + lgamma(nu/2) - (nu/2)*tmp }

tic()
nu <- optimize(Q_nu, interval = c(2 + 1e-16, 1000))$minimum
toc()


# solve optimal nu via bisection method:
# log(nu/2) - digamma(nu/2) = y
optimize_nu <- function(y, tol = 1e-4) {
  if (y < 0) stop("y must be positive!")
  L <- 2 + 1e-6
  U <- 1000

  while (1) {
    mid <- (L + U) / 2
    tmp <- log(mid)-digamma(mid) - y

    if (abs(tmp) < tol) break
    if (tmp > 0) L <- mid else U <- mid
  }

  return(2*mid)
}

tic()
optimize_nu(-1-tmp)
toc()

library(rbenchmark)
benchmark("optimize" = {nu1 <- optimize(Q_nu, interval = c(2 + 1e-16, 1000))$minimum},
          "bisection" = {nu2 <- optimize_nu(-1-tmp)},
          replications = 1000)

Q_nu(nu1)
Q_nu(nu2)
