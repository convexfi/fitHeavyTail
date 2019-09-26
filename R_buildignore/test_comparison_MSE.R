library(mvtnorm)
library(covHeavyTail)
library(ggplot2)
library(reshape2)
#library(latex2exp)

#
# Test of different packages for estimation of mu and R under heavy tails
#
# Comparison with other packages:
#   QRM::fit.mst() for Student-t: fails, not robust
#   MASS::cov.trob(X, nu=nu) for Student-t: it's fast, but requires nu and the performance is not very good
#   MASS::cov.rob(X, method="mcd" or "mve")$cov: it's very slow and the performance not good
#   robustbase::covMcd(X)$cov is slow and not good performance
#   robust::covRob(X)$cov does not perform well
#   covRobust::cov.nnve(X)[c("mu","cov")] does not perform well and it's kind of slow
#   EstimateMoments_MinskerWei(X) is not good...
#   [covHeavyTail or FinCovMat]
#

N <- 20
nu <- 4
mu <- rep(0, N)

set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")


N_realiz <- 20  # multiple realizations for averaging
T_sweep <- round(seq(from = ceiling(2.5*N), to = 5*N, length.out = 20))
MSE_QRM_T <- MSE_Studentt_T <- MSE_Cauchy_T <- MSE_Tyler_T <- MSE_SCM_T <- NULL
time_SCM_T <- NULL
pbar <- txtProgressBar(min = it<-0, max = length(T_sweep), style=3)
for(T in T_sweep) {
  setTxtProgressBar(pbar, it<-it+1)
  MSE_QRM <- MSE_Studentt <- MSE_Cauchy <- MSE_Tyler <- MSE_SCM <- 0
  time_SCM <- 0
  for (ii in 1:N_realiz) {
    # generate data
    X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
    #X <- rmvnorm(n = T, mean = mu, sigma = Sigma)  # Gaussian data

    # Sigma estimation
    time_SCM <- time_SCM + system.time({
      Sigma_SCM <- cov(X)
    })["elapsed"]


    Tyler <- momentsTyler(X)
    Sigma_Tyler <- Tyler$cov

    Cauchy <- momentsCauchy(X)
    Sigma_Cauchy <- Cauchy$cov

    Studentt <- momentsStudentt(X)
    Sigma_Studentt <- Studentt$cov

    Sigma_QRM <- QRM::fit.mst(X, method = "BFGS", nit = 100, tol = 1e-6)$covariance
    #browser()




    # compute MSEs
    MSE_SCM <- MSE_SCM + norm(Sigma_SCM - Sigma, "F")^2
    MSE_Tyler <- MSE_Tyler + norm(Sigma_Tyler - Sigma, "F")^2
    MSE_Cauchy <- MSE_Cauchy + norm(Sigma_Cauchy - Sigma, "F")^2
    MSE_Studentt <- MSE_Studentt + norm(Sigma_Studentt - Sigma, "F")^2
    MSE_QRM  <- MSE_QRM  + norm(Sigma_QRM - Sigma, "F")^2
  }
  # keep track of MSEs
  MSE_SCM_T  <- c(MSE_SCM_T, MSE_SCM / N_realiz)
  MSE_Tyler_T <- c(MSE_Tyler_T, MSE_Tyler / N_realiz)
  MSE_Cauchy_T <- c(MSE_Cauchy_T, MSE_Cauchy / N_realiz)
  MSE_Studentt_T <- c(MSE_Studentt_T, MSE_Studentt / N_realiz)
  MSE_QRM_T <- c(MSE_QRM_T, MSE_QRM / N_realiz)
  # keep track of MSEs
  time_SCM_T <- c(time_SCM_T, time_SCM / N_realiz)
}


# MSE plots
MSE_all_T <- cbind("SCM"        = MSE_SCM_T,
                   "Tyler"      = MSE_Tyler_T,
                   "Cauchy"     = MSE_Cauchy_T,
                   "Student-t"  = MSE_Studentt_T,
                   "QRM"        = MSE_QRM_T)
rownames(MSE_all_T) <- T_sweep
ggplot(melt(MSE_all_T), aes(x = Var1, y = value, col = Var2, shape = Var2)) +
  geom_line() + geom_point() + ylim(0, 200) +
  theme(legend.title = element_blank()) +
  ggtitle(bquote("MSE of covariance matrix estimation for heavy-tailed data (" * N == .(N) * "," ~ nu == .(nu)* ")")) +
  #ggtitle(TeX(sprintf("MSE of covariance matrix estimation for heavy-tailed data (N = %d, $\\nu$ = %.2f)", N, nu))) +
  xlab("T") + ylab("MSE")
# check here is TeX problem has been fixed: https://github.com/stefano-meschiari/latex2exp/pull/21


