set.seed(1234)

# creating constants in 0)
p<-500
pi<-0.5
rho<-0.8
c1<-1
c0<-0

b1<-rep(0,p)
for (i in 1:p) {b1[i]<-1/i}
K<-sqrt(sum(b1^2))/2
b1<-b1/K
sqrt(sum(b1^2))
#b0 is reversed order b1
b0<-rev(b1)

#function generating data
DGP<-function(n){
  # X - design matrix (this is n individuals with known p=500 characteristics)
  library(mvtnorm)
  Sigma <- outer(1:p, 1:p, FUN=function(x, y) rho^(abs(x-y)))
  X <- rmvnorm(n, sigma = Sigma, method = "chol")
  
  # treatment variable W 
  W<-rbinom(n=n,size = 1, prob = pi)
  
  # generate two outcome variables
  eps0<-rnorm(n, 0, 1)
  Y0<-c0+X%*%b0+eps0
  
  eps1<-rnorm(n, 0, 1)
  Y1<-c1+X%*%b1+eps1
  
  #compute tau infeasible
  meanX<-colMeans(X, dims = 1)
  tau_if<-c1-c0+meanX%*%(b1-b0)
  
  # Y
  Y<-W*Y1+(1-W)*Y0
  
  return(list(list(Y,W,X),tau_if))
}

### generate data
#n<-100
#l<-DGP(n)

taufunction<-function(l){

data<-l[[1]]
Y<-data[[1]]
W<-data[[2]]
X<-data[[3]]
tau_if<-as.vector(l[[2]])

## Simple estimator - difference in means
n1<-sum(W)
n0<-sum(1-W)
meanY1<-(1/n1)*W%*%Y
meanY0<-(1/n0)*(1-W)%*%Y
tau_simple<-as.vector(meanY1-meanY0)

# Variance of simple estimator
Var_simple<-(1/n1^2)*W%*%(Y-rep(meanY1,n))^2 + (1/n0^2)*(1-W)%*%(Y-rep(meanY0,n))^2

#confidence interval 90%
upper_simple<-tau_simple+qnorm(p=0.95)*sqrt(Var_simple)
lower_simple<-tau_simple-qnorm(p=0.95)*sqrt(Var_simple)
ci_simple<-c(lower_simple,upper_simple)


#adjusted estimates
#install.packages("devtools")
library(devtools)
#install_github("swager/crossEstimation")
library(crossEstimation)

res <- ate.glmnet(X, Y, W, alpha = 1, conf.level = 0.9, nfolds = 10, method = "separate",
                  lambda.choice = "lambda.min")

return(list(tau_if=tau_if, tau_simple=tau_simple,tau_adj=res$tau,ci_simple=ci_simple,ci_adj=res$conf.int))

}

#T<-taufunction(l)


Indicator<-function(x, min, max)
{ 
  as.numeric(x >= min) * as.numeric(x <= max)
}

Precision<-function(T){
  
SE_1<-(T$tau_if-T$tau_simple)^2
SE_2<-(T$tau_if-T$tau_adj)^2
I1<-Indicator(T$tau_if, T$ci_simple[1],T$ci_simple[2])
I2<-Indicator(T$tau_if, T$ci_adj[1],T$ci_adj[2])

return(list(SE_1=SE_1,SE_2=SE_2,I1=I1,I2=I2))
}

#P<-Precision(T)

#Simulation

Z<-matrix(data=0,nrow = 1000,ncol=4)

n<-100
for (j in 1:1000) {
  
  #l<-DGP(n)
  #T<-taufunction(l)
  P<-Precision(taufunction(DGP(n)))
  Z[j,1]<-P$SE_1
  Z[j,2]<-P$SE_2
  Z[j,3]<-P$I1
  Z[j,4]<-P$I2
  
}

#Results
meanZ<-colMeans(Z, dims = 1, na.rm = TRUE)
meanZ

#SE_1 ; SE_2; I1; I2
#n=50: 0.7672225 0.5061011 0.9530000 0.9663677
#n=100: 0.36396351 0.09567676 0.96800000 0.99900000





