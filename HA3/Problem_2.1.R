library(devtools)
library(optrdd)

set.seed(1234)

# creating constants 
n<-4000
B<-100

tau_true1<-0.2
tau_true2<-0.4
tau_true3<-0

Indicator<-function(x, min, max)
{ as.numeric(x >= min) * as.numeric(x <= max)}


Simulation<-function(){

#Data generation
X = sample(seq(-4, 4, length.out = 41), n, replace = TRUE)
W1<-(X>=0)*1
W2 <- as.numeric(-1 <X & X < 1)
eps<-rnorm(n, 0, sqrt(2))
Y1<-10+0.2*(1+X)*W1+1/(1+exp(2*X))+5*X+eps
Y2<-10+0.2*(1+X)*W2+1/(1+exp(2*X))+5*X+eps



out.1 = optrdd(X=X, Y=Y1, W=W1, max.second.derivative = 0.5, estimation.point = 0)
#print(out.1)
r1<-out.1$tau.hat+out.1$tau.plusminus
l1<-out.1$tau.hat-out.1$tau.plusminus
I1<-Indicator(tau_true1, l1,r1)
se1<-out.1$sampling.se

out.2 = optrdd(X=X, Y=Y2, W=W2, max.second.derivative = 0.5, estimation.point = 1)
#print(out.2) 

r2<-out.2$tau.hat+out.2$tau.plusminus
l2<-out.2$tau.hat-out.2$tau.plusminus
I2<-Indicator(tau_true2, l2,r2)
se2<-out.2$sampling.se

out.3 = optrdd(X=X, Y=Y2, W=W2, max.second.derivative = 0.5, estimation.point = -1)
#print(out.3)

r3<-out.3$tau.hat+out.3$tau.plusminus
l3<-out.3$tau.hat-out.3$tau.plusminus
I3<-Indicator(tau_true3, l3,r3)
se3<-out.3$sampling.se

out.4 = optrdd(X=X, Y=Y2, W=W2, max.second.derivative = 0.5)
r4<-out.4$tau.hat+out.4$tau.plusminus
l4<-out.4$tau.hat-out.4$tau.plusminus
I4<-Indicator(tau_true2, l4,r4)
I5<-Indicator(tau_true3, l4,r4)
se4<-out.4$sampling.se


return(list(I1=I1,se1=se1,I2=I2,se2=se2,I3=I3,se3=se3,I4=I4,se4=se4,I5=I5))
}
Simulation()

#Simulation
Z1<-matrix(data=0,nrow = B,ncol=2)
Z2<-matrix(data=0,nrow = B,ncol=2)
Z3<-matrix(data=0,nrow = B,ncol=2)
Z4<-matrix(data=0,nrow = B,ncol=3)


for (j in 1:B) {
  P<-Simulation()
  Z1[j,1]<-P$I1
  Z1[j,2]<-P$se1
  Z2[j,1]<-P$I2
  Z2[j,2]<-P$se2
  Z3[j,1]<-P$I3
  Z3[j,2]<-P$se3
  Z4[j,1]<-P$I4
  Z4[j,2]<-P$I5
  Z4[j,3]<-P$se4
}
meanZ1<-colMeans(Z1, dims = 1, na.rm = TRUE)
meanZ2<-colMeans(Z2, dims = 1, na.rm = TRUE)
meanZ3<-colMeans(Z3, dims = 1, na.rm = TRUE)
meanZ4<-colMeans(Z4, dims = 1, na.rm = TRUE)


meanZ1
#[1] 0.9300000 0.1774146
meanZ2
#[1] 0.9900000 0.1755997
meanZ3
#[1] 0.9700000 0.1751649
meanZ4
#[1] 0.7900000 0.8300000 0.1324127
