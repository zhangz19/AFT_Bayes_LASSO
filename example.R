

rm(list=ls())
require(survival)
require(flexsurv)
require(glmnet)

set.seed(23)
n <- 500  #or fully 5000 # the size of a cohort
# simulating 20 covariates =
p <- 20 # the number of covariates 
# coefficient vector 
betavec <- c(0.5, 0.5, 0.35, -0.35, rep(0, 16))
z <- matrix(NA, nrow=n, ncol=p)
for(i in 1:p)  z[,i] <- rbinom(n, 1, 0.5)
mu <- z%*%betavec
epsilon.star <- rweibull(n, 0.82, 29.27)
epsilon <- rnorm(n, 0, 1) #log(epsilon.star)/sqrt(1.57)
mytime <- exp( 1 + mu + exp(-0.5*mu^2)*epsilon )
cen.time <- z[,1] + z[,2] + runif(n, 0, 200)
myv <- apply(cbind(mytime, cen.time), 1, min)
mydel <- as.numeric(mytime<cen.time)
#+++++++++ prepare your data. CAUTION: put covariate matrix z at first
dat <- data.frame(z, time=myv, status=mydel)

#************** Cox-LASSO
y <- cbind(time=dat$time, status=dat$status)
x <- as.matrix(dat[,1:p])
cvfit <- cv.glmnet(x, y, family="cox") #first do 10-fold cross-validation to select lambda
m <- glmnet(x, y, family="cox", lambda=cvfit$lambda.min) #plugin the optimal lambda
cat('Cox-LASSO selected: ', paste0('X', which(m$beta!=0)), '\n')

#************** step-genF
forward <- TRUE  #TRUE=forward selection, FALSE=backward elimination
if(forward){ #+++++  forward selection
  fit0 <- flexsurvreg(Surv(time, status) ~ 1, data=dat, dist='lnorm')
  coefs <- as.numeric(coef(fit0))
  init0 <- c(coefs[1],exp(coefs[2]), 0, 1e-6, coefs[-c(1:2)])
  fit <- update(fit0, inits=init0, data=dat, dist='genf')
  myfit <- fit
  len <- length(indc <- 1:p)
  q <- length(inds <- numeric(0))
  aic0 <- BIC(myfit)
  aics <- aic0 - 1
  while(q<=p){
    aic0 <- BIC(myfit)
    aics <- numeric(len)
    fits <- list(len)
    for(i in 1:len){
      inds0 <- c(inds, indc[i])
      myform0 <- as.formula(paste("Surv(time, status) ~ ", paste(names(dat)[inds0], collapse= "+")))
      fit0 <- flexsurvreg(myform0, data=dat, dist='lnorm',method="Nelder-Mead")
      coefs <- as.numeric(coef(fit0))
      init0 <- c(coefs[1],exp(coefs[2]),0, 1e-6, coefs[-c(1:2)])
      fits[[i]] <- update(fit0, inits=init0, dist='genf',method="Nelder-Mead")
      aics[i] <- BIC(fits[[i]])
    }
    optind <- which(aics == min(aics))
    if(min(aics) < aic0){
      myfit <- fits[[optind]]
      # cat('X',indc[optind],' was included.\n',sep='')
      q <- length( inds <- c(inds, indc[optind]) )
      len <- length(indc <- indc[-optind])
    }else q <- p+1 #stop
  }
}else{ #----- backward selection
  myform0 <-  as.formula(paste("Surv(time, status) ~ ", paste(names(dat)[1:p], collapse= "+")))
  fit0 <- flexsurvreg(myform0, data=dat, dist='lnorm',method = "Nelder-Mead")
  coefs <- as.numeric(coef(fit0))
  init0 <- c(coefs[1],exp(coefs[2]),0,1e-6, coefs[-c(1:2)])
  fit <- update(fit0, inits=init0, data=dat, dist='genf')
  myfit <- fit
  q <- length(inds <- 1:p)
  aic0 <- BIC(myfit)
  aics <- aic0-1
  while(q>0){
    aic0 <- BIC(myfit)
    aics <- numeric(q)
    fits <- list(q)
    for(i in 1:q){
      inds0 <- inds[-i]
      myform0 <-  as.formula(paste("Surv(time, status) ~ ", paste(names(dat)[inds0], collapse= "+")))
      fit0 <- flexsurvreg(myform0, data=dat, dist='lnorm',method = "Nelder-Mead")
      coefs <- as.numeric(coef(fit0))
      init0 <- c(coefs[1],exp(coefs[2]),0, 1e-6, coefs[-c(1:2)])
      fits[[i]] <- update(fit0, inits=init0, dist='genf')
      aics[i] <- BIC(fits[[i]])
    }
    optind <- which(aics == min(aics))
    if(min(aics) < aic0){
      myfit <- fits[[optind]]
      # cat('X',inds[optind],' was excluded.\n',sep='')
      q <- length( inds <- inds[-optind] )
    }else q <- 0 #stop
  }
}

cat('step-genF selected: ', paste0('X', sort(inds)), '\n')

# p-val
# suppose myfit is a flexsurvreg object (we used dist='genf') 
a <- myfit$res
# in a, column 1 is estimate of beta, column 4 is its associated standard error
# p-value is upper tail probability of the approximated standard Normal density
pval <- ifelse(is.null(nrow(a)), 1-pnorm(abs(a[1]/a[4])), 1-pnorm(abs(a[,1]/a[,4])))
print(cbind(a,pval))

# not run

