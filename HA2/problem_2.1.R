library(reshape2)
library(Synth)
library(glmnet)

set.seed(1234)

# creating constants in 0)
n<-50
T<-35
B<-100

# X and X_L matrices
X <- matrix(10*rnorm(n*T),nrow = n)
svd_x <- svd(X)
X_L <- (svd_x$u)%*%
  (diag(ifelse(abs(svd_x$d)>quantile(abs(svd_x$d),0.7),svd_x$d,0)))%*%
  t(svd_x$v)
X_L <- X_L[sample(1:n,n),]

errorfunction<-function(){
    
  # data generation
  eps<-rnorm(n, 0, sqrt(3))
  
  F0<-matrix(0,nrow = n, ncol=T)
  for (t in 1:T) {F0[,t]<-0.5*t+sin(t)}
  F1<-matrix(0,nrow = n, ncol=T)
  for (t in 1:T) { F1[,t]<-0.1*max(t-20,0)+0.2*(max(t-20,0))^2 }
  W<-matrix(0,nrow = n, ncol=T)
  for (t in 21:T) { W[1,t]<-1 }
  
  EPS <- matrix(rnorm(n*T, mean = 0, sd = sqrt(3)),nrow = n)
  
  Y<-X_L+F0+F1*W+EPS
  
  
  
  #Transform data to panel like type
  
  Y_df <- data.frame(Y) 
  colnames(Y_df)<-c(1:T)
  Y_df$id <- 1:nrow(Y_df) 
  Y_df <- melt(Y_df, id.vars = 'id')
  colnames(Y_df)<-c("i","t","Y")
  
  
  
  
  #This is special command for Synth command, that generates Z0 and Z1
  #we specify predictors as t variable, because the function need something for predictors
  #then will replace X0 and X1 with Z0 and Z1 respectively
  
  dataprep.out <- dataprep(
    foo = Y_df,
    dependent = "Y",
    predictors = "t",
    unit.variable = "i",
    time.variable = "t",
    treatment.identifier = 1,
    controls.identifier = c(2:n),
    time.predictors.prior = 1:20,
    time.optimize.ssr = 1:20,
    time.plot = 1:T)
  
  dataprep.out$X0<-dataprep.out$Z0
  dataprep.out$X1<-dataprep.out$Z1
  
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS",custom.v = rep(1/20,20))
  
  #### THIS IS TREATMENT EFFECT OVER TIME for SC!
  tau_sc <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w)
  
  ## below are some tables and graphs (not asked in the task) about the resulting estimates
  synth.tables <- synth.tab(dataprep.res = dataprep.out,
                            synth.res = synth.out)
  path.plot(synth.res = synth.out, dataprep.res = dataprep.out,
            Ylab = "Y", Xlab = "period",
            Ylim = c(-25, 40), Legend = c("Treated",
                                        "synthetic control"), Legend.position = "bottomright")
  gaps.plot(synth.res = synth.out, dataprep.res = dataprep.out,
            Ylab = "gap in Y", Xlab = "period",
            Ylim = c(-40, 40), Main = NA)
  #computing the error between estimated and actual treatment
  tau_actual<-as.data.frame(t(F1[1,]))
  t(tau_sc)
  colnames(tau_actual)<-c(1:T)
  error_sc<-abs(tau_actual-tau_sc)/tau_actual
  error_sc<-error_sc[,c(21:T)]
  
  
  ## Second estimator
  
  #y is pretreatment (i.e. periods 1-20) outcome for the unit 1
  #x is pretreatment (i.e. periods 1-20) outcomes for the unit 2-35
  TY<-t(Y)
  y<-TY[c(1:35),1]
  x<-TY[c(1:35),c(2:n)]
  
  fit <- cv.glmnet(x=x, y=y, alpha = 0.5, lambda =(1:100)/1000, family = "gaussian", foldid = c(1,rep(2, 19), rep(3, 15)), keep=TRUE, lambda.min.ratio=0.0001)
  diff <- (fit$fit.preval - y%*%t(rep(1, 100)))
  mse <- colSums(diff[21:35, 1:100]^2)
  plot(fit$lambda, mse)
  bestlambda<-fit$lambda[which.min(mse[fit$lambda > 0.02])]
  
  newx<-TY[c(21:T),c(2:n)]
  tau_en<-predict(fit, newx = newx, s = bestlambda)
  
  
  #computing the error between estimated and actual treatment
  tau_actual_short<-as.data.frame(t(F1[1,c(21:T)]))
  tau_en<-t(tau_en)
  colnames(tau_actual_short)<-c(21:T)
  tau_en<-as.data.frame(tau_en)
  colnames(tau_en)<-c(21:T)
  
  error_en<-abs(tau_actual_short-tau_en)/tau_actual_short
  
  return(list(error_en=error_en,error_sc=error_sc))
}


EN<-as.data.frame(matrix(data=0,nrow = B,ncol=15))
SC<-as.data.frame(matrix(data=0,nrow = B,ncol=15))
for (j in 1:B) {
  P<-errorfunction()
  EN[j,]<-P$error_en[1,]
  SC[j,]<-P$error_sc[1,]
}

error_en<-colMeans(EN, dims = 1, na.rm = TRUE)
error_sc<-colMeans(SC, dims = 1, na.rm = TRUE)

plot(error_sc,  xlab="Period", ylab="Error ", pch=19, lty=1)
lines(error_sc,  xlab="Period", ylab="Error ", pch=19, lty=1)
lines(error_en, col="red",lty=2)
points(error_en, col="red", pch="*")
legend("topright",legend=c("Sinthetic Control","Elastic net"), col=c("black","red"),
       pch=c(".","*"),lty=c(1,2), ncol=1)


