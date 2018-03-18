library(mvtnorm)
library(zeallot)
library(crossEstimation)

set.seed(1234)


p <- 500
pi <- 0.5
rho <- 0.8
c_1 <- 1
c_0 <- 0

beta_1 <-  1/seq(1, p)
K = norm(beta_1, type = "2")/2
beta_1 <- beta_1/K
beta_0 <- rev(beta_1)

generate_data <- function(){
  Sigma <- outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  
  W <- rbinom(n = n, size = 1, prob = pi)
  
  eps_0 <- rnorm(n)
  Y0 <- c_0 + X %*% beta_0 + eps_0
  eps_1 <- rnorm(n)
  Y1 <- c_1 + X %*% beta_1 + eps_1
  
  meanX <- colMeans(X)
  tau_if <- c_1 - c_0 + meanX %*% (beta_1 - beta_0)
  
  Y <- W * Y1 + (1-W) * Y0
  
  return(list(data=list(Y=Y, W=W, X=X), tau_if=tau_if))
}

estimate_effect <- function(data){
  c(dataset, tau_if) %<-% data
  c(Y, W, X) %<-% dataset
  
  ## Simple estimator - difference in means
  n_1 <- sum(W)
  n_0 <- sum(1-W)
  meanY1 <- (1/n_1) * sum(W*Y)
  meanY0 <- (1/n_0) * sum((1-W) * Y)
  tau_simple <- meanY1 - meanY0
  
  # Variance of simple estimator
  Var_simple <- (1/n_1^2) * sum(W * (Y - meanY1)^2) +
      (1/n_0^2) * sum((1-W) * (Y - meanY0)^2)
  
  #confidence interval 90%
  upper_simple <- tau_simple + qnorm(p=1 - ((1 - 0.90)/2)) * sqrt(Var_simple)
  lower_simple <- tau_simple - qnorm(p=1 - ((1 - 0.90)/2)) * sqrt(Var_simple)
  ci_simple <- c(lower_simple, upper_simple)
  
  
  #adjusted estimates
  res <- ate.glmnet(X, Y, W, alpha = 1, conf.level = 0.9, nfolds = 10, method = "separate",
                    lambda.choice = "lambda.min")
  
  return(list(tau_if = tau_if,
              tau_simple = tau_simple,
              tau_adj = res$tau,
              ci_simple = ci_simple,
              ci_adj = res$conf.int))
}



evatuate_estimates <- function(estimates){
  c(tau_if, tau_simple, tau_adj, ci_simple, ci_adj) %<-% estimates
  
  se_1 <- (tau_if - tau_simple)^2
  se_2 <- (tau_if - tau_adj)^2
  I_1 <- tau_if >= ci_simple[1] && tau_if <= ci_simple[2]
  I_2 <- tau_if >= ci_adj[1] && tau_if <= ci_adj[2]

return(c(se_1, se_2, I_1, I_2))
}

simulation_instance <- function(place_holder) {
  data <- generate_data()
  params <- estimate_effect(data)
  return(evatuate_estimates(params))
}

simulation <- function(iterations) {
  Z <- matrix(0, nrow = iterations, ncol=4)
  Z <- pbapply(Z, 1, simulation_instance)
  return(t(Z))
}

assess_simulation <- function() {
  result <- numeric()
  iterations <- c(50, 100, 300, 500, 1000)
  
  for (it in iterations) {
    op <- pboptions(type = "timer")
    Z <- simulation(it)
    pboptions(op)
    result <- rbind(result, colMeans(Z))
  }
  return(result)
}

assess_simulation()
