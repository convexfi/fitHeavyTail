library(fitHeavyTail)
library(tictoc)

set.seed(234)

N <- 100
T <- 4*N
factors <- 5
beta <- matrix(rnorm(N*factors, 10, 1), N, factors)
psi <- rexp(N, 1/10)
Sigma_scale <- beta %*% t(beta) + diag(psi)  # scale matrix (not covariance)
mu <- rnorm(N, 0, 1)
nu <- 7

X <- mvtnorm::rmvt(n = T, sigma = Sigma_scale, df = nu, delta = mu)

# a function for adding NAs
addNA <- function(Y, missing_ratio) {
  n <- nrow(Y)
  p <- ncol(Y)
  n_missing <- round(n*missing_ratio)
  missing_mask <- sample(n, n_missing)

  Y_incomp <- Y
  for (i in 1:n_missing) {
    mask_NA <- sample(c(TRUE, FALSE), p, replace = TRUE, prob = c(0.1, 0.9))
    Y_incomp[missing_mask[i], mask_NA] <- NA
  }
  return(Y_incomp)
}
X_wNA <- addNA(X, 0.1)



tic()
fit_old <- momentsStudentt(X)
toc()

# Naming:
# fit_mvt <- momentsStudentt
# fit_mvtFA <- covTFA
# fit_mvskewt


tic()
fit_nom <- fit_mvt(X, ftol = 1e6, return_iterates = TRUE)  # use return_convergence? return_iterations? Or debug?
toc()

method <- "ECM"
fit_old <- fit_mvt(X, ptol = 0, ftol = 0, return_iterates = TRUE)
fit_nom <- fit_mvt(X, ptol = 0, ftol = 0, return_iterates = TRUE)

fit_wFA <- fit_mvt(X, factors = factors, ftol = 1e6, return_iterates = TRUE)
fit_wFA_wNA <- fit_mvt(X_wNA, factors = factors, ftol = 1e6, return_iterates = TRUE)
# fit_wFA <- covTFA(X, factors = factors, ptol = ptol, ftol = ftol, return_iterates = TRUE)  # check better name for factors?
# fit_wFA_wNA <- covTFA(X_wNA, factors = factors, ptol = ptol, ftol = ftol, return_iterates = TRUE)

fit_old$nu
fit_nom$nu

norm(fit_old$mu - mu, "2") / norm(mu, "2")
norm(fit_nom$mu - mu, "2") / norm(mu, "2")
norm(fit_old$mu - fit_nom$mu, "2") / norm(fit_nom$mu, "2")
norm(fit_old$cov * (fit_old$nu - 2) / fit_old$nu - fit_nom$Sigma, "F") / norm(fit_nom$Sigma, "F")


plot(sapply(fit_nom$iterations_record, function(x) x$nu))
sapply(fit_wFA$iterations_record, function(x) x$nu)
sapply(fit_wFA_wNA$iterations_record, function(x) x$nu)


# estimate nu from subspace of X ---------------------------
library(mvtnorm)
library(fitHeavyTail)
library(ggplot2)

N <- 10
T <- 1.1*N
nu <- 6
mu <- rep(0, N)

# generate data
set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")

X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
apply(X, 2, function(x) fit_mvt(cbind(x))$nu)
median(sapply(as.list(1:10), function(x) fit_mvt(X[, sample(N, 2)], ptol = 1e-3)$nu))

res_mvt_noreg <- fit_mvt(X, return_iterates = TRUE)
res_mvt_wreg  <- fit_mvt(X, nu_regcoef = 2e-2, return_iterates = TRUE)

p_mvt_noreg <- fitHeavyTail:::plotConvergence(res_mvt_noreg)
do.call(gridExtra::grid.arrange, c(p_mvt_noreg, ncol = 1))

p_mvt_wreg <- fitHeavyTail:::plotConvergence(res_mvt_wreg)
do.call(gridExtra::grid.arrange, c(p_mvt_wreg, ncol = 1))


# test for univariate case ---------------------------
N <- 1
T <- 10
nu <- 10
mu <- rep(0, N)

# generate data
set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")

X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
res_mvt_noreg <- fit_mvt(X, return_iterates = TRUE)



# verify if the regularization term works --------------------------
N <- 20
nu <- 20
mu <- rep(0, N)

# set.seed(357)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma
qplot(eigen(Sigma)$values, geom = "histogram", xlab = "eigenvalues", fill = I("cyan"), col = I("black"),
      main = "Histogram of eigenvalues of true covariance matrix")

N_realiz <- 5  # multiple realizations for averaging
T_sweep <- round(seq(from = ceiling(1.5*N), to = 2.5*N, length.out = 8))
MSE_Studentt_T  <- MSE_Studentt_nu3_T <- MSE_Studentt_nu4_T <- MSE_Studentt_nu5_T <- MSE_Studentt_nu6_T <- MSE_Studentt_nu7_T <- MSE_Studentt_normal_T  <- MSE_SCM_T <- NULL
time_Studentt_T <- time_Studentt_nu_T <- time_SCM_T <- NULL

pbar <- txtProgressBar(min = it<-0, max = length(T_sweep), style=3)
for(T in T_sweep) {
  setTxtProgressBar(pbar, it<-it+1)
  MSE_Studentt  <- MSE_Studentt_nu3 <- MSE_Studentt_nu4 <- MSE_Studentt_nu5 <- MSE_Studentt_nu6 <- MSE_Studentt_nu7 <- MSE_Studentt_normal <- MSE_SCM  <- 0
  time_Studentt <- time_Studentt_nu <- time_SCM <- 0
  for (ii in 1:N_realiz) {
    # generate data
    X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)  # heavy-tailed data
    #X <- rmvnorm(n = T, mean = mu, sigma = Sigma)  # Gaussian data


    # Sigma estimation
    time_SCM <- time_SCM + system.time({
      Sigma_SCM <- cov(X)
    })["elapsed"]


    Sigma_Studentt_nu3 <- fit_mvt(X, nu = 3)$cov
    Sigma_Studentt_nu4 <- fit_mvt(X, nu = 4)$cov
    Sigma_Studentt_nu5 <- fit_mvt(X, nu = 5)$cov
    Sigma_Studentt_nu6 <- fit_mvt(X, nu = 6)$cov
    Sigma_Studentt_nu7 <- fit_mvt(X, nu = 7)$cov


    Sigma_Studentt_normal <- fit_mvt(X)$cov

    time_Studentt <- time_Studentt + system.time({
      Sigma_Studentt <- fit_mvt(X, nu_target = 6, nu_regcoef = 5e-2)$cov
      # res <- fit_mvt(X, ftol = 1e6, return_iterates = TRUE)
      # plot(sapply(res$iterations_record, `[[`, "nu"))
      # plot(sapply(res$iterations_record, `[[`, "log_likelihood"))
    })["elapsed"]


    # compute MSEs
    MSE_SCM <- MSE_SCM + norm(Sigma_SCM - Sigma, "F")^2
    MSE_Studentt_nu3 <- MSE_Studentt_nu3 + norm(Sigma_Studentt_nu3 - Sigma, "F")^2
    MSE_Studentt_nu4 <- MSE_Studentt_nu3 + norm(Sigma_Studentt_nu4 - Sigma, "F")^2
    MSE_Studentt_nu5 <- MSE_Studentt_nu3 + norm(Sigma_Studentt_nu5 - Sigma, "F")^2
    MSE_Studentt_nu6 <- MSE_Studentt_nu6 + norm(Sigma_Studentt_nu6 - Sigma, "F")^2
    MSE_Studentt_nu7 <- MSE_Studentt_nu3 + norm(Sigma_Studentt_nu7 - Sigma, "F")^2
    MSE_Studentt_normal <- MSE_Studentt_normal + norm(Sigma_Studentt_normal - Sigma, "F")^2
    MSE_Studentt <- MSE_Studentt + norm(Sigma_Studentt - Sigma, "F")^2
 }
  # keep track of MSEs
  MSE_SCM_T  <- c(MSE_SCM_T, MSE_SCM / N_realiz)
  MSE_Studentt_nu3_T <- c(MSE_Studentt_nu3_T, MSE_Studentt_nu3 / N_realiz)
  MSE_Studentt_nu4_T <- c(MSE_Studentt_nu4_T, MSE_Studentt_nu4 / N_realiz)
  MSE_Studentt_nu5_T <- c(MSE_Studentt_nu5_T, MSE_Studentt_nu5 / N_realiz)
  MSE_Studentt_nu6_T <- c(MSE_Studentt_nu6_T, MSE_Studentt_nu6 / N_realiz)
  MSE_Studentt_nu7_T <- c(MSE_Studentt_nu7_T, MSE_Studentt_nu7 / N_realiz)
  MSE_Studentt_normal_T <- c(MSE_Studentt_normal_T, MSE_Studentt_normal / N_realiz)
  MSE_Studentt_T <- c(MSE_Studentt_T, MSE_Studentt / N_realiz)

}

# MSE plots
MSE_all_T <- cbind("SCM"                   = MSE_SCM_T,
                   "Student-t (nu=3)"      = MSE_Studentt_nu3_T,
                   "Student-t (nu=4)"      = MSE_Studentt_nu4_T,
                   "Student-t (nu=5)"      = MSE_Studentt_nu5_T,
                   "Student-t (nu=6)"      = MSE_Studentt_nu6_T,
                   "Student-t (nu=7)"      = MSE_Studentt_nu7_T,
                   "Student-t normal"      = MSE_Studentt_normal_T,
                   "Student-t"             = MSE_Studentt_T)
rownames(MSE_all_T) <- T_sweep
ggplot(melt(MSE_all_T), aes(x = Var1, y = value, col = Var2, shape = Var2)) +
  geom_line() + geom_point() + coord_cartesian(ylim = c(0, 1000)) +
  theme(legend.title = element_blank()) +
  ggtitle(bquote("MSE of covariance matrix estimation for heavy-tailed data (" * N == .(N) * "," ~ nu == .(nu)* ")")) +
  xlab("T") + ylab("MSE")


# check the logic for treating nu ----------------------------------------
rm(list = ls())
library(mvtnorm)
library(fitHeavyTail)
library(ggplot2)
library(reshape2)

N <- 20
T <- 40
nu <- 5
mu <- rep(0, N)

set.seed(123)
U <- t(rmvnorm(n = round(1.1*N), sigma = 0.1*diag(N)))
Sigma <- U %*% t(U) + diag(N)
Sigma_scale <- (nu-2)/nu * Sigma

X <- rmvt(n = T, delta = mu, sigma = Sigma_scale, df = nu)

# original MLE without any tricks on nu
fit_mvt(X)$nu

# fix nu to a give vlaue, say 6
fit_mvt(X, nu = 6)$nu

# fix nu to a give vlaue, but decided via kurtosis
fit_mvt(X, nu = "kurtosis")$nu

# regularize nu to a given target
fit_mvt(X, nu_target = 6, nu_regcoef = 1)$nu
fit_mvt(X, nu_target = 6, nu_regcoef = 1e10)$nu

fit_mvt(X, nu_regcoef = 1)$nu
fit_mvt(X, nu_regcoef = 1e10)$nu

# sanity check for simultaneously passing nu and nu_target
fit_mvt(X, nu = 3, nu_target = 6, nu_regcoef = 1)$nu
fit_mvt(X, nu = 3, nu_target = 6, nu_regcoef = 1e10)$nu

