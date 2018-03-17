library(mvtnorm)

set.seed(1234)

n <- 100
p <- 500
rho <- 0.5
tau <- 0.5

beta_Y <-  1/(seq(1, p)^(2))
beta_Y <- beta_Y/norm(beta_Y, type = "2")

beta_W <- 1/(seq(1, p)^(1/2))
beta_W <- beta_W/norm(beta_W, type = "2")


generate_data <- function() {
  Sigma <- outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  
  theta <- X %*% beta_W + rnorm(n)
  W <- rbinom(n, 1, 1/(1+exp(theta)))
  
  epsilon <- rnorm(n)
  Y <- tau*W + X %*% beta_Y + epsilon
  
  return(list(Y, X, W))
}
