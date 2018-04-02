library('mvtnorm')
library('glmnet')
library('stats')
library('rpart')
library('zeallot')
library('pbapply')

set.seed(1234)

n_E <- 500
n_O <- 500
B <- 200

p_set <- c(5, 10, 20, 30, 50, 70, 100, 130, 170, 200)

generate_parameters <- function(p) {
  alpha_S <- rnorm(p, mean = 0, sd = 1/sqrt(p))
  gamma_S <- alpha_S
  return(list(alpha_S=alpha_S, gamma_S=gamma_S))
}

generate_dataset <- function(alpha_S, gamma_S, n_O, n_E, p) {
  
  e_X <- function(a, s){
    output <- as.numeric(exp(a %*% s) / (1 + exp(a %*% s)))
    return(output)
  }
  
  Ip <- diag(x = 1, p, p)
  S_O <- rmvnorm(n=n_O, sigma = Ip, method = "chol")
  S_E <- rmvnorm(n=n_E, sigma = Ip, method = "chol")
  
  e_WE <- e_X(alpha_S, t(S_E)) 
  e_WO <- e_X(alpha_S, t(S_O))
  
  e_YE <- e_X(gamma_S, t(S_E))
  e_YO <- e_X(gamma_S, t(S_O))
  
  w_E <- rbinom(n_E, 1, e_WE)
  w_O <- rbinom(n_E, 1, e_WO)
  Y_E <- rbinom(n_E, 1, e_YE)
  Y_O <- rbinom(n_O, 1, e_YO)
  return(list(S_O = S_O,
              S_E = S_E,
              w_E = w_E,
              w_O = w_O,
              Y_E = Y_E,
              Y_O = Y_O))
}

estimate <- function(data) {
  c(S_O, S_E, w_E, w_O, Y_E, Y_O) %<-% data
  h_function.fit <- cv.glmnet(S_O, Y_O, alpha=0.5, family="gaussian", nfolds=10, lambda.min.ratio=0.01)
  Y_E.hat <- as.numeric(predict(h_function.fit, S_E, type="response"))
  
  g_1.hat <- w_E / sum(w_E)
  g_0.hat <- (1-w_E) / sum((1-w_E))
  
  tau_surr.hat <- sum(Y_E.hat*g_1.hat) - sum(Y_E.hat*g_0.hat)
  tau_exp.hat <- sum(Y_E*g_1.hat) - sum(Y_E*g_0.hat)
  
  SE <- (tau_surr.hat - tau_exp.hat)^2
  return(SE)
}

simulation_instance <- function(params) {
  c(alpha_S, gamma_S, n_O, n_E, p) %<-% params
  data <- generate_dataset(alpha_S, gamma_S, n_O, n_E, p)
  SE <- estimate(data)
  return(SE)
}

simulation <- function(iterations, params) {
  Z <- rep(list(params), iterations)
  Z <- pbsapply(Z, simulation_instance)
  return(Z)
}

assess_simulation <- function() {
  res_Z <- matrix(nrow=B)
  for (p in p_set) {
    print(paste0("p = ", p))
    c(alpha_S, gamma_S) %<-% generate_parameters(p)
    params <- list(alpha_S=alpha_S, gamma_S=gamma_S, n_O=n_O, n_E=n_E, p=p)
    op <- pboptions(type = "timer")
    Z <- simulation(B, params)
    pboptions(op)
    res_Z <- cbind(res_Z, Z)
  }
  res_Z <- res_Z[,2:dim(res_Z)[2]]
  colnames(res_Z) <- p_set
  return(res_Z)
}

SE <- assess_simulation()
avgSE <- colMeans(SE)
plot(avgSE, type="l", xaxt = "n", ylab = "squarred error", xlab = "number of surrogates")
axis(1, at=1:10, labels=p_set)
