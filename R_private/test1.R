library(mvtnorm)
library(devtools)
library(rbokeh)
library(covHeavyTail)

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

N <- 40
#T <- 200
set.seed(357)
nv <- 4
mu <- rep(0,N)
#u <- rnorm(N, sd=0.1)
#R_cov <- u %o% u + diag(N)
U <- t(rmvnorm(n=round(0.7*N), sigma=0.1*diag(N)))
R_cov <- U %*% t(U) + diag(N)
R_scale <- (nv-2)/nv*R_cov   # R_cov=nv/(nv-2)*R_scale
hist(eigen(R_cov)$values, breaks=200, prob=TRUE, col="cyan", main="Histogram of covariance matrix")


T_sweep <- seq(from=50, by=40, to=2000)
error_Cauchy <- error_Tyler <- error_Student_t <- error_NSCM <- error_SCM <- NULL
pbar <- txtProgressBar(min = it<-0, max = length(T_sweep), style=3)
for(T in T_sweep) {
  setTxtProgressBar(pbar, it<-it+1)

  #generate data
  #X <- rmvnorm(n=T, mean=mu, sigma=R_cov)  #Gaussian data
  X <- mvtnorm::rmvt(n=T, delta=mu, sigma=R_scale, df=nv)  #heavy-tailed data


  # SCM estimation
  R_SCM <- cov(X)
  #X_ <- scale(X, center=TRUE, scale=FALSE)
  #R_SCM <- (1/(nrow(X)-1)) * t(X_) %*% X_


  #NSCM
  row_norms <- sqrt(apply(X^2, 1, "sum"))
  X_norm <- sweep(X, 1, FUN="/", STATS=row_norms)
  R_NSCM <- (N/T) * t(X_norm) %*% X_norm
  #R_NSCM <- N*cov(X_norm)


  #Heavy-tailed estimators
  Cauchy <- momentsCauchy(X, verbose=FALSE)
  Student_t <- momentsStudentt(X, verbose=TRUE)
  Tyler <- momentsTyler(X, verbose=FALSE)



  #compute errors
  error_SCM <- c(error_SCM, norm(R_SCM - R_cov, "F"))
  error_NSCM <- c(error_NSCM, norm(R_NSCM - R_cov, "F"))
  error_Student_t <- c(error_Student_t, norm(Student_t$cov - R_cov, "F"))
  error_Tyler <- c(error_Tyler, norm(Tyler$cov - R_cov, "F"))
  error_Cauchy <- c(error_Cauchy, norm(Cauchy$cov - R_cov, "F"))
  # cat(sprintf( "[T=%d, N=%d] Errors:  SCM: %g;  NSCM: %g;  Student-t: %g;  K-tau: %g;  R-Tyler: %g;  R-Cauchy: %g;  R-MinskerWei: %g\n", T, N,
  #              last(error_SCM),
  #              last(error_NSCM),
  #              last(error_Student_t),
  #              last(error_Tyler),
  #              last(error_Cauchy) ))
}

#plots
print( figure(width=800, title=sprintf("Estimation error of covariance matrix (N=%d)",N), xlab="T", ylab="error", legend_location = "top_right") %>%
         ly_lines(T_sweep, error_SCM, color="blue", legend="SCM")  %>%
         ly_lines(T_sweep, error_Tyler, color="orange", width=3, legend="Median-Tyler") %>%
         ly_lines(T_sweep, error_Cauchy, color="purple", width=1, legend="Cauchy") %>%
         ly_lines(T_sweep, error_Student_t, color="red", width=2, legend="Student-t") )



