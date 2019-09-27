library(mvtnorm)
library(covHeavyTail)
library(ggplot2)
library(reshape2)

#
# Test of different packages for estimation of mu and R under heavy tails
#
# Comparison with other packages:
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
T_sweep <- round(seq(from = ceiling(1.5*N), to = 5*N, length.out = 12))
MSE_MASS_cov.trov_T <- MSE_QRM_T <- MSE_Studentt_T <- MSE_Cauchy_T <- MSE_Tyler_T <- MSE_SCM_T <- NULL
time_MASS_cov.trov_T <- time_QRM_T <- time_Studentt_T <- time_Cauchy_T <- time_Tyler_T <- time_SCM_T <- NULL
pbar <- txtProgressBar(min = it<-0, max = length(T_sweep), style=3)
for(T in T_sweep) {
  setTxtProgressBar(pbar, it<-it+1)
  MSE_MASS_cov.trov <- MSE_QRM <- MSE_Studentt <- MSE_Cauchy <- MSE_Tyler <- MSE_SCM <- 0
  time_MASS_cov.trov <- time_QRM <- time_Studentt <- time_Cauchy <- time_Tyler <- time_SCM <- 0
  for (ii in 1:N_realiz) {
    # generate data
    X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
    #X <- rmvnorm(n = T, mean = mu, sigma = Sigma)  # Gaussian data


    # Sigma estimation
    time_SCM <- time_SCM + system.time({
      Sigma_SCM <- cov(X)
    })["elapsed"]

    time_Tyler <- time_Tyler + system.time({
      Sigma_Tyler <- momentsTyler(X)$cov
    })["elapsed"]

    time_Cauchy <- time_Cauchy + system.time({
      Sigma_Cauchy <- momentsCauchy(X)$cov
    })["elapsed"]

    time_Studentt <- time_Studentt + system.time({
      Sigma_Studentt <- momentsStudentt(X)$cov
    })["elapsed"]

    time_QRM <- time_QRM + system.time({
      Sigma_QRM <- tryCatch(
        as.matrix(QRM::fit.mst(X, method = "BFGS", nit = 50, tol = 1e-6)$covariance),
        warning = function(w) return(NA),
        error   = function(e) return(NA))
    })["elapsed"]

    time_MASS_cov.trov <- time_Studentt + system.time({
      Sigma_MASS_cov.trov <- MASS::cov.trob(X, nu = 4)$cov
    })["elapsed"]





    # compute MSEs
    MSE_SCM <- MSE_SCM + norm(Sigma_SCM - Sigma, "F")^2
    MSE_Tyler <- MSE_Tyler + norm(Sigma_Tyler - Sigma, "F")^2
    MSE_Cauchy <- MSE_Cauchy + norm(Sigma_Cauchy - Sigma, "F")^2
    MSE_Studentt <- MSE_Studentt + norm(Sigma_Studentt - Sigma, "F")^2
    MSE_QRM  <- MSE_QRM  + norm(Sigma_QRM - Sigma, "F")^2
    MSE_MASS_cov.trov  <- MSE_MASS_cov.trov  + norm(Sigma_MASS_cov.trov - Sigma, "F")^2
  }
  # keep track of MSEs
  MSE_SCM_T  <- c(MSE_SCM_T, MSE_SCM / N_realiz)
  MSE_Tyler_T <- c(MSE_Tyler_T, MSE_Tyler / N_realiz)
  MSE_Cauchy_T <- c(MSE_Cauchy_T, MSE_Cauchy / N_realiz)
  MSE_Studentt_T <- c(MSE_Studentt_T, MSE_Studentt / N_realiz)
  MSE_QRM_T <- c(MSE_QRM_T, MSE_QRM / N_realiz)
  MSE_MASS_cov.trov_T <- c(MSE_MASS_cov.trov_T, MSE_MASS_cov.trov / N_realiz)
  # keep track of MSEs
  time_SCM_T <- c(time_SCM_T, time_SCM / N_realiz)
  time_Tyler_T <- c(time_Tyler_T, time_Tyler / N_realiz)
  time_Cauchy_T <- c(time_Cauchy_T, time_Cauchy / N_realiz)
  time_Studentt_T <- c(time_Studentt_T, time_Studentt / N_realiz)
  time_QRM_T <- c(time_QRM_T, time_QRM / N_realiz)
  time_MASS_cov.trov_T <- c(time_MASS_cov.trov_T, time_MASS_cov.trov / N_realiz)
}


# MSE plots
MSE_all_T <- cbind("SCM"            = MSE_SCM_T,
                   "Tyler"          = MSE_Tyler_T,
                   "Cauchy"         = MSE_Cauchy_T,
                   "Student-t"      = MSE_Studentt_T,
                   "QRM"            = MSE_QRM_T,
                   "MASS::cov.trov" = MSE_MASS_cov.trov_T)
rownames(MSE_all_T) <- T_sweep
ggplot(melt(MSE_all_T), aes(x = Var1, y = value, col = Var2, shape = Var2)) +
  geom_line() + geom_point() + coord_cartesian(ylim = c(0, 500)) +
  theme(legend.title = element_blank()) +
  ggtitle(bquote("MSE of covariance matrix estimation for heavy-tailed data (" * N == .(N) * "," ~ nu == .(nu)* ")")) +
  xlab("T") + ylab("MSE")
# ggtitle(latex2exp::TeX(sprintf("MSE of covariance matrix estimation for heavy-tailed data (N = %d, $\\nu$ = %.2f)", N, nu))) +
# check here is TeX problem has been fixed: https://github.com/stefano-meschiari/latex2exp/pull/21


# time plots
time_all_T <- cbind("SCM"            = time_SCM_T,
                    "Tyler"          = time_Tyler_T,
                    "Cauchy"         = time_Cauchy_T,
                    "Student-t"      = time_Studentt_T,
                    "QRM"            = time_QRM_T,
                    "MASS::cov.trov" = time_MASS_cov.trov_T)
rownames(time_all_T) <- T_sweep
time_all_T[is.na(MSE_all_T)] <- NA
ggplot(melt(time_all_T), aes(x = Var1, y = value, col = Var2, shape = Var2)) +
  geom_line() + geom_point() +
  theme(legend.title = element_blank()) +
  ggtitle(bquote("Computational cost of the computation of covariance matrices for heavy-tailed data (" * N == .(N) * "," ~ nu == .(nu)* ")")) +
  xlab("T") + ylab("cpu time")
