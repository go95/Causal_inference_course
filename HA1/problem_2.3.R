#Random forest

set.seed(1234)

# creating constants in 0)
n<-1000
p<-10
pi<-0.5
X.test = matrix(0, 101, p)
X.test[,1] = seq(-2, 2, length.out = 101)

#function generating data
DGP<-function(n){
  # X - design matrix 
  library(mvtnorm)
  Sigma <- diag(p)
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  
  # treatment variable W 
  W<-rbinom(n=n,size = 1, prob = pi)
  
  # generate  outcome variable
  eps<-rnorm(n, 0, 1)
  Y<-pmax(0,X[,1])*W+X[,2]+pmax(0, X[,2])+eps
  
  return(list(Y=Y,W=W,X=X))
}
#data<-DGP(n)


taufunction<-function(data){
library(grf)
tau.forest <- causal_forest(data$X, data$Y, data$W, num.trees = 4000)
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)

tau_hat<-tau.hat$predictions
CI<- list(L=tau.hat$predictions - qnorm(0.95) * sigma.hat, U=tau.hat$predictions + qnorm(0.95) * sigma.hat)

#true effects
tau_true<-pmax(0,X.test[,1])

return(list(tau_true=tau_true,tau_hat=tau_hat,CI=CI))
}

#R<-taufunction(data)

Indicator<-function(x, min, max)
{ 
  as.numeric(x >= min) * as.numeric(x <= max)
}

Errorfunction<-function(R){
  SE=(R$tau_true-R$tau_hat)^2
  I<-matrix(Indicator(R$tau_true, R$CI$L,R$CI$U), nrow = 101, ncol = 1)
  return(cbind(SE,I))
}

#Errorfunction(R)


#Simulation
k=300

Z<-matrix(data=0,nrow = 101,ncol=k)
I<-matrix(data=0,nrow = 101,ncol=k)
for (j in 1:k) {
  
  P<-Errorfunction(taufunction(DGP(n)))
  Z[,j]<-P[,1]
  I[,j]<-P[,2]
}

meanZ<-rowMeans(Z, dims = 1, na.rm = TRUE)
meanI<-rowMeans(I, dims = 1, na.rm = TRUE)

plot(X.test[,1],y =  meanZ,  xlab="x", ylab="Error ", pch=19, main = "p=10")

plot(X.test[,1],y =  meanI,  xlab="x", ylab="Coverage ", pch=19, main = "p=10")


