

library('mvtnorm')
library('grf')
library('devtools')
library('crossEstimation')
library('balanceHD')


# =================================================================
# Problem 2.1 ------------- Experimental Simulation ---------------
# =================================================================

# 0
# Global Parameters

set.seed(1234)

n <- 100
p <- 500
pi <- 0.5
rho <- 0.8
c1 <- 1
c0 <- 0

b1 <- seq(1, p, by = 1)^(-1)
K <- sqrt(sum(b1^2))/2
b1 <- b1/K
b0 <- rev(b1)


# ------- [1] ------- 

Generator <- function(n){
  Sigma <- outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  W <- rbinom(n=n, size = 1, prob = pi)
  e0 <- rnorm(n, mean=0, sd=1)
  Y0 <- c0 + X%*%b0 + e0
  e1 <- rnorm(n, 0, 1)
  Y1 <- c1 + X%*%b1 + e1
  avg.X <- colMeans(X, dims = 1)
  tau.if <- c1 - c0 + avg.X %*% (b1 - b0)
  Y <- W*Y1 + (1-W)*Y0
  data <- list(Y=Y, W=W, X=X)
  output <- list(Data=data, Tau.IF=tau.if)
  return(output)
}


# ------- [2] ------- 

data <- Generator(n)

F <- function(data){
  Y <- data$Data$Y
  X <- data$Data$X
  W <- data$Data$W
  y0 <- Y*W
  Y0 <- y0[y0!=0]
  y1 <- Y*(1-W)
  Y1 <- y1[y1!=0]
  avgY0 <- mean(Y0)
  avgY1 <- mean(Y1)
  tau.simple <- avgY0 - avgY1
  varY0 <- var(Y0)
  varY1 <- var(Y1)
  V = varY0 + varY1
  alpha <- 0.9
  q <- qnorm(1-(1-alpha)/2)
  CI.simple_low <- tau.simple - q * sqrt(V)
  CI.simple_up <- tau.simple + q * sqrt(V)
  CI.simple <- list(CI_low = CI.simple_low, CI_up = CI.simple_up)
  res <- ate.glmnet(X, Y, W, alpha = 1, nfolds = 10, method = 'separate',
                    lambda.choice = 'lambda.min')
  tau.adj <- res$tau
  CI.adj <- res$conf.int
  CI.adj_low <- CI.adj[1]
  CI.adj_up <- CI.adj[2]
  CI.adj <- list(CI_low = CI.adj_low, CI_up = CI.adj_up)
  output <- list(Tau.Simple = tau.simple, CI.Simple = CI.simple, Tau.Adj = tau.adj, CI.Adj = CI.adj)
  return(output)
}


# ------- [3] ------- 

taus <- F(data)

SqErr <- function(taus){
  tau.if <- data$Tau.IF
  tau.simple <- taus$Tau.Simple
  tau.adj <- taus$Tau.Adj
  SE1 <- (tau.if - tau.simple)^2
  SE2 <- (tau.if - tau.adj)^2
  ci.simple <- taus$CI.Simple
  ci.adj <- taus$CI.Adj
  I1 <- tau.if >= ci.simple$CI_low & tau.if < ci.simple$CI_up
  I2 <- tau.if >= ci.adj$CI_low & tau.if < ci.adj$CI_up
  output <- c(SE1, SE2, I1, I2)
  return(output)
}


# ------- [4] ------- 









# ------- [5 (Bonus) ] -------



# =================================================================
# Problem 2.2 ------------ Ovservational Simulation ---------------
# =================================================================

# 0
# Global Parameters

set.seed(1234)

n <- 100
p <- 500
tau <- 0.5
rho <- 0.8

BY <- seq(1, p, by = 1)^(-2)
BW <- seq(1, p, by = 1)^(-1/2)
cy <- sqrt(sum(BY^2))
cw <- sqrt(sum(BW^2))
BY <- BY/cy
BW <- BW/cw


# ------- [1] ------- 

Generator <- function(n){
  Sigma <- outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  theta <- X %*% bw + rnorm(n)
  W <- rbinom(n, 1, 1/(1+exp(theta)))
  e <- rnorm(n, mean=0, sd=1)
  Y <- tau*W + X%*%by + e
  output <- list(Y, W, X)
  names(output) <- c('Y', 'W', 'X')
  return(output)
}


# ------- [2] ------- 

data <- Generator(n)

F <- function(data){
  Y <- data$Data$Y
  X <- data$Data$X
  W <- data$Data$W
  y0 <- Y*W
  Y0 <- y0[y0!=0]
  y1 <- Y*(1-W)
  Y1 <- y1[y1!=0]
  avgY0 <- mean(Y0)
  avgY1 <- mean(Y1)
  tau.simple <- avgY0 - avgY1
  varY0 <- var(Y0)
  varY1 <- var(Y1)
  V = varY0 + varY1
  tau.DL <- rlassoEffect(X, Y, W, method = 'double selection')
  tau.DB <- residualBalance.ate(X, Y, W, target.pop = 1,
                                fit.method = 'elnet', alpha = 0.9, zeta = 0.5)
  output <- list(Tau.DL = tau.DL, Tau.DB = tau.DB, Tau.Simple = tau.simple)
  return(output)
}


# ------- [3] ------- 

taus <- F(data)

SqErr <- function(taus){
  tau.dl <- taus$Tau.DL
  tau.db <- taus$Tau.DB
  tau.simple <- taus$Tau.Simple
  err.dl <- (tau - tau.dl)^2
  err.db <- (tau - tau.db)^2
  err.simple <- (tau - tau.simple)^2
  output <- list(SE.DL = err.dl, SE.DB = err.db, SE.Simple = err.simple)
}


# ------- [4] ------- 




# ------- [5 (Bonus) ] -------





# =================================================================
# Problem 2.3 --------- Heterogeneous Treatment Effects -----------
# =================================================================

# 0
# Global Parameters

set.seed(1234)

n <- 1000
p <- 10
pi <- 0.5

X.test <- matrix(0, 101, p)
X.test[,1] <- seq(-2, 2, length.out = 101)


# ------- [1] ------- 

Generator <- function(n){
  X <- rmvnorm(n, sigma = diag(p), method = "chol")
  W <- rbinom(n = n, size = 1, prob = pi)
  errs <- rnorm(n, mean = 0, sd = 1)
  Y <- pmax(X[,1], 0) * W + X[,2] + pmax(X[,2], 0) + errs
  return(list(Y=Y, W=W, X=X))
}


# ------- [2] ------- 

data <- Generator(n)


F <- function(data){
  tau.forest <- causal_forest(data$X, data$Y, data$W, num.trees = 4000)
  tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
  t.hat <- tau.hat$predictions
  sigma.hat <- sqrt(tau.hat$variance.estimates)
  alpha <- 0.9
  q <- qnorm(1-(1-alpha)/2)
  CI_low <- t.hat - q * sigma.hat
  CI_up <- t.hat + q * sigma.hat
  CI <- list(CI_low = CI_low, CI_up = CI_up)
  tau.true <- pmax(X.test[,1], 0)
  return(list(Tau.True = tau.true, Tau.Hat = t.hat, CI = CI))
}

# ------- [3] ------- 

taus <- F(data)

SqErr <- function(taus){
  SE <- (taus$Tau.True - taus$Tau.Hat) ^ 2
  I <- taus$Tau.True >= taus$CI$CI_low & taus$Tau.True <= taus$CI$CI_up
  M <- as.matrix(cbind(SE, I))
  colnames(M) <- c("SEs", "Coverage")
  return(M)
}


# ------- [4] ------- 




# ------- [5 (Bonus) ] -------





