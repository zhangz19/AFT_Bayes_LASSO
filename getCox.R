

# Cox-LASSO with bootstrap for simulation data
# no bootstrap: 500 simulations are completed with CPUtime 29.613 minutes!

rm(list=ls())

sumPro <- FALSE

ptm <- proc.time()[3]

#id <- as.numeric(commandArgs(TRUE))   # id=1:nbin

# # some interpretation
#https://www.medcalc.org/manual/cox_proportional_hazards.php

## results for proportional hazard model LASSO
library(R.matlab)
library(survival)
library(glmnet)
library(c060)
library(peperr)
library(methods) # key for cv.glmnet, see https://github.com/dhimmel/elevcan/issues/1

# simulation design in the paper
load('t0')
len <- length(t0)
Z0 <- matrix(c(0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1), ncol=4)
Z0[1,] <- c(1,1,1,0)
x0 <- matrix(0, 4, 20);
x0[,1:4] <- Z0
ncas <- nrow(Z0)

dirs <- dir(patter='simudat*')
nsim <- length(dirs)

# now split nsim jobs into bins
nbin <- 500
starts <- seq(1,nsim,by=nsim/nbin)
ends <- starts+nsim/nbin-1

mydirs <- './'; #'save/'
dirs <- dir(mydirs, patter='simudat*')

nboot <- 500
mysBeta <- mymBeta <- array(0,c(nsim,20))
mat <- array(NA, c(len, ncas, nsim))

mystart <- 1
#for(k in mystart:myend){
for(k in 1:nsim){
  cat(k,''); if(!k%%20) cat('\n')
  #print(k)
  sdat <- readMat(paste(mydirs,'simudat',k,'.mat',sep=''))
  
  nn <- nrow(sdat$Z)
  
  # get estimates
  x <- as.matrix(sdat$Z)
  y <- cbind(time=as.numeric(sdat$V), status=as.integer(sdat$Delta))
  cvfit <- cv.glmnet(x,y, family="cox")
  fit <- glmnet(x=x,y=y, family="cox", lambda=cvfit$lambda.min)
  
  # this might be a bug in that package
  fit$response <- y
  fit$linear.predictor <- x%*%fit$beta
  
  a <- predictProb(object=fit, response=y, times=t0, x=x0,complexity=cvfit$lambda.min)
  
  #link = x0%*%fit$beta, response = exp(x0%*%fit$beta)
  mysBeta[k-mystart+1,which(fit$beta!=0)] <- 1
  mymBeta[k-mystart+1,] <- as.numeric(fit$beta)
  mat[,,k-mystart+1] <- t(a)
  fit0 <- fit
  beta0 <- as.numeric(fit0$beta)
  indsig <- which(beta0!=0)
  
  # # get bootstrap confidence intervals
  #for(k_boot in 1:nboot){
  #cat(k_boot,''); if(!k_boot%%20) cat('\n')
  #
  # inds <- sample(1:nn, nn, replace=T)# sample (Y,X), we treat X as random
  #
  # x <- as.matrix(sdat$Z[inds,])
  # y <- cbind(time=as.numeric(sdat$V[inds]), status=as.integer(sdat$Delta[inds]))
  # 
  # cvfit <- cv.glmnet(x,y, family="cox")
  # fit <- glmnet(x=x,y=y, family="cox", lambda=cvfit$lambda.min)
  # 
  # #fit <- coxph(Surv(y[,1], y[,2])~x[,indsig])
  # 
  # #y0 <- y; y0[,1] <- round(y0[,1])
  # 
  # # this might be a bug in that package
  # fit1$response <- y
  # fit1$linear.predictor <- x[,indsig]%*%fit$coefficients
  # 
  # a <- predictProb(object=fit1, response=y, times=t0, x=x0,complexity=cvfit$lambda.min)
  # 
  # # # this is not predicted survival probability
  # # a <- predict(fit, newx=x0, type='response') # link'response'
  # 
  # #link = x0%*%fit$beta, response = exp(x0%*%fit$beta)
  # mysBeta[k-mystart+1,which(fit$beta!=0),1+k_boot] <- 1
  # mat[,,k-mystart+1,1+k_boot] <- t(a)
  #}
  
}

cat('\n')

cputime <- as.numeric(proc.time()[3]-ptm)/60
cat(nsim, 'simulations are completed with CPUtime', round(cputime,3), 'minutes!','\n')

save(file='Csp_s1',mat=mat, mysBeta=mysBeta,mymBeta=mymBeta)


if(sumPro){
  nbin <- 500
  Z0 <- matrix(c(0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1), ncol=4)
  Z0[1,] <- c(1,1,1,0)
  ncas <- nrow(Z0)
  x0 <- matrix(0, 4, 20);
  x0[,1:4] <- Z0
  dirs <- dir(patter='simudat*')
  nsim <- length(dirs)
  nboot <- 500
  load('t0')
  len <- length(t0)
  mysBetaAll <- array(0,c(nsim,20, 1+nboot))
  matAll <- array(NA, c(len, ncas, nsim, 1+nboot))
  starts <- seq(1,nsim,by=nsim/nbin)
  ends <- starts+nsim/nbin-1  
  for(i in 1:nbin){
    load(paste('outCox',i,sep=''))
    mystart <- starts[i]; myends <- ends[i]
    mysBetaAll[mystart:myends,,] <- mysBeta[1:(myends-mystart+1),,]
    matAll[,,mystart:myends,] <- mat[,,1:(myends-mystart+1),]
  }  
  mat <- matAll
  mysBeta <- mysBetaAll
  mat <- mat[,,,1]
  mysBeta <- mysBeta[,,1]
  save(file='Csp_s2',mat=mat, mysBeta=mysBeta)
}

