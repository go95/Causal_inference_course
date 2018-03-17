library(mvtnorm)
library(zeallot)
library(hdm)
library(balanceHD)

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

estimate_effect <- function(data) {
  c(Y, X, W) %<-% data
  
  Eff <- rlassoEffect(X, Y, W, method = "double selection")
  tau_dl <- coefficients(Eff)
  
  tau_db <- residualBalance.ate(X, Y, W, target.pop = 1,
                                fit.method = "elnet", alpha = 0.9, zeta = 0.5)
  
  tau_simple <- mean(X[!W]) - mean(X[W])
  
  return(c(tau_dl, tau_db, tau_simple))
}


squared_error <- function(params) {
  return((params - tau)^2)
}

simulation_instance <- function(place_holder) {
  data <- generate_data()
  params <- estimate_effect(data)
  return(squared_error(params))
}

simulation <- function(iterations) {
  Z <- matrix(0, nrow=iterations, ncol=3)
  Z <- apply(Z, 1, simulation_instance)
  return(t(Z))
}

assess_simulation <- function() {
  Z <- simulation(400)
  return(colMeans(Z))
}

main()
