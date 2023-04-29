library(microbenchmark)
library(ggplot2)


nu <- 4
N <- 200
T <- 600
T_sweep  <- seq(from = 30, to = 90, by = 10)
nu_sweep <- seq(from = 4, to = 12, by = 2)

# covmat from real market data
X <- diff(log(SP500_2015to2020$stocks))[-1, ]
mu_big            <- colMeans(X)
Sigma_cov_big     <- cov(X)
Sigma_scatter_big <- (nu-2)/nu * Sigma_cov_big

# extract parameters with correct dimensions
set.seed(42)
idx <- sample(ncol(X), N)
mu            <- mu_big[idx]
Sigma_cov     <- Sigma_cov_big[idx,idx]
Sigma_scatter <- (nu-2)/nu * Sigma_cov

X <- rmvt(n = T, delta = mu, sigma = Sigma_scatter, df = nu)


# sanity check
r2_a <- diag( X %*% inv(Sigma_scatter) %*% t(X) )
r2_b <- rowSums(X * (X %*% solve(Sigma_scatter)))
r2_c <- rowSums(X * t(solve(Sigma_scatter, t(X))))
all.equal(r2_a, r2_b)
all.equal(r2_a, r2_c)


res <- microbenchmark(
  "naive"          = diag( X %*% inv(Sigma_scatter) %*% t(X) ),
  "inverse matrix" = rowSums(X * (X %*% solve(Sigma_scatter))),
  "linear system"  = rowSums(X * t(solve(Sigma_scatter, t(X)))),
  setup = {X <- rmvt(n = T, delta = mu, sigma = Sigma_scatter, df = nu)},
  times = 1000 #, check = "equal"
)


print(res)
plot(res)



p <- ggplot(res, aes(x = time/1000, y = expr, color = expr)) +  # Plot performance comparison
  geom_violin() +
  geom_boxplot(width = 0.1) #+ scale_x_continuous(trans = 'log2')

print(p)
