library(mvtnorm)
library(grf)
library(zeallot)
library(pbapply)

set.seed(1234)

n <- 1000
p <- 20 #this is just for compatability with bonus task
pi <- 0.5
X.test <- matrix(0, 101, p)
X.test[,1] <- seq(-2, 2, length.out = 101)

generate_data <- function() { 
  Sigma <- diag(p)
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  W <- rbinom(n=n,size = 1, prob = pi)
  
  # generate  outcome variable
  eps <- rnorm(n)
  Y <- pmax(0, X[,1]) * W + X[,2] + pmax(0, X[,2]) + eps
  
  return(list(Y = Y, W = W, X = X))
}


estimate_effect <- function(data, p){

  tau.forest <- causal_forest(data$X[,1:p], data$Y, data$W, num.trees = 4000)
  tau.hat <- predict(tau.forest, X.test[,1:p], estimate.variance = TRUE)
  sigma.hat <- sqrt(tau.hat$variance.estimates)
  
  tau_hat <- tau.hat$predictions
  CI <- c(tau.hat$predictions - qnorm(0.95) * sigma.hat, tau.hat$predictions + qnorm(0.95) * sigma.hat)
  
  tau_true <- pmax(0,X.test[,1])
  
  return(list(tau_true = tau_true, tau_hat = tau_hat, CI = CI))
}


squared_error <- function(estimates) {
  se = (estimates$tau_true - estimates$tau_hat)^2
  I <- estimates$tau_true >= estimates$CI[1] && estimates$tau_true <= estimates$CI[1]
  return(c(se, I))
}



simulation_instance <- function(place_holder, p) {
  data <- generate_data()
  params <- estimate_effect(data, p)
  return(squared_error(params))
}

simulation <- function(iterations, p) {
  Z <- matrix(0, nrow = iterations, ncol=3)
  Z <- pbapply(Z, 1, simulation_instance, p=p)
  return(t(Z))
}


assess_simulation <- function() {
  for p in c(1, 3, 20) {
    result <- numeric()
    op <- pboptions(type = "timer")
    Z <- simulation(300, p)
    pboptions(op)
    means <- colMeans(Z)
    plot(X.test[,1], result)
    result <- rbind(result, means)
  }
  return(result)
}

result <- assess_simulation()
