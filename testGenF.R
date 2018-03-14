


rm(list=ls())

runChunk <- 1

if(runChunk == 1){
  ##----------- chunk 1: run in HPC for different repl
  options(warn=-1)  #disable warnings
  
  # argument repl indicates the bin number
  repl <- as.numeric(commandArgs(TRUE))
  
  myfilename <- paste('out',repl,sep='')
  
  # generalized F
  library(R.matlab)
  library(flexsurv)
  library(survival)
  
  M <- 5e2
  mysBeta_gf <- matrix(0,M,20)
  mydirs_gf = './'
  
  forward <- FALSE
  mycriterion <- "BIC"
  myalpha <- - 0.05
  mydist <- c('lnormal','gengamma','genf')
  useLnormal <- TRUE
  
  select_criterion <- function(fit, opt=c("AIC","BIC","p-value"), 
                               forward=TRUE, useLnormal=FALSE, verbose=FALSE){ 
    if(opt == "AIC") return(AIC(fit))
    if(opt == "BIC") return(BIC(fit))
    if(opt == "p-value"){
      # methods for p-value based backwards selection
      # fit must be a flexsurvreg fit with ldist=genf
      nullvalue <- 1
      if(forward) nullvalue <- 1
      t0 <- 4; if(useLnormal) t0 <- 2
      a <- fit$res[-c(1:t0), ]   #4 for genf
      if(length(a)>0){
        if(is.null(nrow(a))) pvals <- 1-pnorm(abs(a[1]/a[4]))
        else pvals <- 1-pnorm(abs(a[,1]/a[,4]))
        ntrial <- length(mytol <- c(1e-16, 1e-22, 1e-30, 1e-40, 1e-60, 1e-65, 1e-70, 1e-80))
        mytrial <- 1
        while(any(is.na(pvals))){
          if(verbose) cat('\nrefit trial = ', mytrial,'...', sep='')
          fit <- update(fit, inits=fit$res[,1], control=list(fnscale=2500, reltol=mytol[mytrial]))  #
          a <- fit$res[-c(1:t0), ]   #4 for genf
          if(is.null(nrow(a))) pvals <- 1-pnorm(abs(a[1]/a[4]))
          else pvals <- 1-pnorm(abs(a[,1]/a[,4]))
          mytrial <- mytrial + 1
          if(any(is.na(pvals)) && mytrial==(ntrial+1)){
            print("Hessian failed")  # stop(" Hessian failed")
            pvals <- nullvalue #forget about this case
            break
          }
        }
        optind <- as.integer(which(pvals == max(pvals)))
        if(forward) optind <- length(pvals) # the last one: new variable
        if(!forward && all(pvals==nullvalue)){ 
          optind <- length(pvals) # remove the last one
          pvals[-optind] <- 0.5 
          print("Hessian failed, sacrifice the last variable")
        }
        return(list(optind=optind, pval= - pvals[optind]))
      }
      else return(list(optind=0, pval= -nullvalue))
    } 
  }
  
  len <- length(t0 <- seq(0,60,by=.2))  # need to be consistent!
  # simulation in the paper
  Z0 <- matrix(c(0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1), ncol=4)
  Z0[1,] <- c(1,1,1,0)
  x0 <- matrix(0, 4, 20);
  x0[,1:4] <- Z0
  ncas <- nrow(Z0)
  mat <- array(NA, c(len, ncas, M))
  
  ptm <- proc.time()[3]
  count <- 0
  
  for(k in repl){  #1:M
    cat('\n\nreplication =',k,'\n')
    print(repl)
    if(!k%%20) cat('\n')
    sdat <- readMat(paste(mydirs_gf,'simudat',k,'.mat',sep=''))
    p <- ncol(tmp <- as.data.frame(sdat$Z))
    names(tmp) <- paste('X',1:p,sep='')
    tmp$time <- as.numeric(sdat$V)
    tmp$status <- as.integer(sdat$Delta) #verify
    
    if(forward){###
      # forward selection
      fit0 <- flexsurvreg(Surv(time, status) ~ 1, data=tmp, dist='lnorm', control=list(fnscale = 2500))
      coefs <- as.numeric(coef(fit0))
      init0 <- c(coefs[1],exp(coefs[2]),0,1e-6, coefs[-c(1:2)])
      fit <- fit0
      if(!useLnormal) fit <- update(fit0, inits=init0, data=tmp, dist='genf')
      # fit <- flexsurvreg(Surv(time, status) ~ 1, data=tmp, dist='genf', control=list(fnscale = 2500))
      
      myfit <- fit
      len <- length(indc <- 1:p)
      q <- length(inds <- numeric(0))
      aic0 <- select_criterion(myfit, mycriterion, forward, useLnormal)
      if(mycriterion == "p-value"){
        aics <- - aic0$pval
        aic0 <- - myalpha;
      } 
      aics <- aic0-1
      
      while(min(aics) < aic0 && q<=p){
        cat('q = ',q,':  ',sep='')           
        if(mycriterion != "p-value"){
          aic0 <- select_criterion(myfit, mycriterion, forward, useLnormal)
          cat(round(aic0),': ',sep='')
        }
        aics <- rep(1, q)
        fits <- list()
        
        for(i in 1:len){
          cat('X',indc[i],' ',sep='')
          inds0 <- c(inds, indc[i])
          myform0 <-  as.formula(paste("Surv(time, status) ~ ", paste(names(tmp)[inds0], collapse= "+")))
          fit0 <- flexsurvreg(myform0, data=tmp, dist='lnorm',method = "Nelder-Mead", control=list(fnscale=2500, reltol=1e-10))
          coefs <- as.numeric(coef(fit0))
          init0 <- c(coefs[1],exp(coefs[2]),0, 1e-6, coefs[-c(1:2)])
          fits[[i]] <- fit0
          if(!useLnormal) fits[[i]] <- update(fit0, inits=init0, dist='genf')
          aic1 <- select_criterion(fits[[i]], mycriterion, forward, useLnormal)
          if(mycriterion == "p-value") aic1 <- - aic1$pval
          aics[i] <- aic1
          count <- count+1
          cat(round(aics[i],4),', ',sep='')
          if(mycriterion == "p-value"){if(aic1 < 1e-8) break}
        }
        cat('\n')
        optind <- which(aics == min(aics))[1]
        if(min(aics) < aic0){
          myfit <- fits[[optind]]
          cat('X',indc[optind],' was chosen\n\n',sep='')
          q <- length( inds <- c(inds, indc[optind]) )
          len <- length(indc <- indc[-optind])
        }
      }
      
    }###
    
    if(!forward){###
      # backward selection
      myform0 <-  as.formula(paste("Surv(time, status) ~ ", paste(names(tmp)[1:p], collapse= "+")))
      fit0 <- flexsurvreg(myform0, data=tmp, dist='lnorm',method = "Nelder-Mead")
      coefs <- as.numeric(coef(fit0))
      init0 <- c(coefs[1],exp(coefs[2]),0,1e-6, coefs[-c(1:2)])
      fit <- fit0
      if(!useLnormal) fit <- update(fit0, inits=init0, data=tmp, dist='genf')
      
      myfit <- fit
      aic0 <- select_criterion(myfit, mycriterion, forward, useLnormal)
      q <- length(inds <- 1:p)
      if(mycriterion == "p-value"){
        vec <- aic0$optind # only one variable to remove 
        aics <- aic0$pval
        aic0 <- myalpha;
      } 
      else aics <- aic0-1
      
      while(min(aics) < aic0 && q>0){
        cat('q = ',q,' ',sep='')
        if(mycriterion != "p-value") {aic0 <- select_criterion(myfit, mycriterion, forward, useLnormal); cat(round(aic0),': ',sep=''); vec <- 1:q}
        if(mycriterion == "p-value") cat(round(-aics,4),': ',sep='')
        aics <- rep(NA, q)
        fits <- list()
        
        for(i in vec){
          cat('X',inds[i],' ',sep='')
          inds0 <- inds[-i]
          myform0 <-  as.formula(paste("Surv(time, status) ~ ", 
                                       paste(names(tmp)[inds0], collapse= "+")))
          fit0 <- flexsurvreg(myform0, data=tmp, dist='lnorm',method = "Nelder-Mead")
          fits[[i]] <- fit0
          if(!useLnormal){
            coefs <- as.numeric(coef(fit0))
            init0 <- c(coefs[1],exp(coefs[2]),0, 1e-10, coefs[-c(1:2)])
            fits[[i]] <- update(fit0, inits=init0, dist='genf')
          }
          aic1 <- select_criterion(fits[[i]], mycriterion, forward, useLnormal)
          if(mycriterion == "p-value"){
            optind <- aic1$optind
            aic1 <- aic1$pval
            fits1 <- fits[[i]]
          }
          aics[i] <- aic1
          count <- count+1
          if(mycriterion != "p-value") cat(round(aics[i],4),', ',sep='')
          if(mycriterion == "p-value") cat(round(-aics[i],4),', ',sep='')
        }
        cat('\n')
        aics <- as.numeric(na.omit(aics))
        
        if(mycriterion != "p-value"){
          optind <- which(aics == min(aics))
          fits1 <- fits[[optind]]
          vec <- optind
        }
        if(min(aics) < aic0){
          myfit <- fits1
          cat('X',inds[vec],' was excluded\n\n',sep='')
          q <- length( inds <- inds[-vec] )
          vec <- optind
        }
      }
      
    }###
    
    if(mycriterion == "p-value"){
      cat('final model:\n')
      a <- myfit$res
      if(is.null(nrow(a))) pval <- 1-pnorm(abs(a[1]/a[4]))
      else pval <- 1-pnorm(abs(a[,1]/a[,4]))
      print(cbind(a,pval))
    }
    
    for(j in 1:4)  
      mat[,j,repl] <- summary(object=myfit, X=x0[j,inds], type='survival', t=t0)[[1]]$est
    
    mysBeta_gf[k,inds] <- 1
  }
  
  
  cputime <- as.numeric(proc.time()[3]-ptm)
  cputime <- cputime/60
  
  cat(count, 'search steps done with CPUtime', round(cputime,3), 'minutes: completed!','\n')
  save(file=myfilename,'mysBeta_gf','mat')
  
}

if(runChunk==2){
  ##----------- chunk 2: combine the results run in parallel
  len1 <- length(t0 <- seq(0,60,by=.2))  # need to be consistent!
  
  Z0 <- matrix(c(0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1), ncol=4)
  Z0[1,] <- c(1,1,1,0)
  x0 <- matrix(0, 4, 20);
  x0[,1:4] <- Z0
  ncas <- nrow(Z0)
  mycase <- 12
  p <- 20
  M <- 500
  mydir0 <- './'
  len <- length(tmp <- dir(path=mydir0,pattern='out'))
  print(min(M,len))
  mysBeta_gfAll <- matrix(0,min(M,len),p)
  matAll <- array(NA, c(len1, ncas, M)) 
  for(i in 1:min(M,len)){
    load(paste(mydir0,tmp[i],sep='')) 
    k <- as.integer(substr(tmp[i], start=4, stop=10L))
    mysBeta_gfAll[i,] <- mysBeta_gf[k,]
    matAll[,,i] <- mat[,,k]
  }
  
  mysBeta_gf <- na.omit(mysBeta_gfAll)
  mat <- na.omit(matAll)
  print(dim(mysBeta_gf))
  sb3 <- mysBeta_gf
  
  falsePostive3 <- rowSums(sb3[,5:ncol(sb3)])
  falseNegatives3 <- 4-rowSums(sb3[,1:4])
  bw3.dat <- as.data.frame(c(falsePostive3, falseNegatives3)); names(bw3.dat) <- 'y'
  bw3.dat$grp <- as.factor(c(rep('False Positives',length(falsePostive3)),
                             rep('False Negatives',length(falseNegatives3))))
  mean(bw3.dat$y[bw3.dat$grp == 'False Positives']/16) 
  mean(bw3.dat$y[bw3.dat$grp == 'False Negatives']/4)
  
  bw.dat <- bw3.dat
  
  save(file=paste('c',mycase, sep=''),'mysBeta_gf','bw.dat')
  
  # combine data in different folders: work*
  mycase <- 10
  p <- 20
  M <- 500
  mysBeta_gfAll <- matrix(NA,M,p)
  mydir0 <- './' 
  len <- length(tmp <- dir(path=mydir0, pattern='work'))
  
  len1 <- length(t0 <- seq(0,60,by=.2))  # need to be consistent!
  Z0 <- matrix(c(0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1), ncol=4)
  Z0[1,] <- c(1,1,1,0)
  x0 <- matrix(0, 4, 20);
  x0[,1:4] <- Z0
  ncas <- nrow(Z0)
  matAll <- array(NA, c(len1, ncas, M)) 
  
  m <- 1
  for(i in 1:len){
    nam <- paste(mydir0,tmp[i],sep='')
    len2 <- length(tmp2 <- dir(path=nam,pattern='out'))
    for(j in 1:len2){
      load(paste(nam,tmp2[j],sep='/'))
      k <- as.integer(substr(tmp2[j], start=4, stop=10L))
      mysBeta_gfAll[m,] <- mysBeta_gf[k,]
      matAll[,,m] <- mat[,,k]
      m <- m+1
      if(m>M) break
    }
  }
  
  mysBeta_gf <- na.omit(mysBeta_gfAll)
  mat <- na.omit(matAll)
  print(dim(mysBeta_gf))
  
  
  sb3 <- mysBeta_gf
  
  falsePostive3 <- rowSums(sb3[,5:ncol(sb3)])
  falseNegatives3 <- 4-rowSums(sb3[,1:4])
  bw3.dat <- as.data.frame(c(falsePostive3, falseNegatives3)); names(bw3.dat) <- 'y'
  bw3.dat$grp <- as.factor(c(rep('False Positives',length(falsePostive3)),
                             rep('False Negatives',length(falseNegatives3))))
  mean(bw3.dat$y[bw3.dat$grp == 'False Positives']/16) 
  mean(bw3.dat$y[bw3.dat$grp == 'False Negatives']/4)
  
  bw.dat <- bw3.dat
  save(file=paste0('c',mycase),'mysBeta_gf','bw.dat','mat')
  
  # get COX-LASSO
  mycase <- 13
  load('Csp_s2')
  sb3 <- mysBeta
  
  falsePostive3 <- rowSums(sb3[,5:ncol(sb3)])
  falseNegatives3 <- 4-rowSums(sb3[,1:4])
  bw3.dat <- as.data.frame(c(falsePostive3, falseNegatives3)); names(bw3.dat) <- 'y'
  bw3.dat$grp <- as.factor(c(rep('False Positives',length(falsePostive3)),
                             rep('False Negatives',length(falseNegatives3))))
  
  mean(bw3.dat$y[bw3.dat$grp == 'False Positives']/16) 
  mean(bw3.dat$y[bw3.dat$grp == 'False Negatives']/4)
  
  mysBeta_gf <- mysBeta
  bw.dat <- bw3.dat
  
  save(file=paste('./genF/c',mycase, sep=''),'mysBeta_gf','bw.dat')

  # get the Bayesian methods provieded in the paper
  mycase <- 14
  tmp <- readMat('tab_s1.mat')
  trues <- c(0.5,0.5,0.35,-0.35)
  sb <- tmp$sBeta
  falsePostive <- rowSums(sb[,5:ncol(sb)])
  falseNegatives <- 4-rowSums(sb[,1:4])
  bw.dat <- as.data.frame(c(falsePostive, falseNegatives)); names(bw.dat) <- 'y'
  bw.dat$grp <- as.factor(c(rep('False Positives',length(falsePostive)),
                            rep('False Negatives',length(falseNegatives))))
  
  mean(bw.dat$y[bw.dat$grp == 'False Positives']/16) 
  mean(bw.dat$y[bw.dat$grp == 'False Negatives']/4)
  
  mysBeta_gf <- sb
  
  par(las=1, mar=c(1.4,2.2,.4,0)+.4,mgp=c(1.8,.3,0), tck=-0.01, cex.axis=1, cex.lab=1, cex.main=1) 
  boxplot(y~grp, data=bw.dat, ylab='Frequency')
  
  save(file=paste('./genF/s2/c',mycase, sep=''),'mysBeta_gf','bw.dat')
  
  
  postscript(file=paste0('./bar20_s2.ps'),pointsize=12,
             width=5.5,height=4.5,horizontal=F,paper='special')  
  express <- function(char.expressions){
    return(parse(text=paste(char.expressions,collapse=";")))
  } 
  par(las=1, mar=c(1,1.4,.4,0)+.4,mgp=c(1.5,.4,0), tck=-0.01, 
      cex.axis=1, cex.lab=1, cex.main=1, xaxt='n')
  barplot(sb/nrow(sb), col='gray80', border='gray80') 
  text((1:ncol(sb))*1.21-.6,rep(-0.03,ncol(sb)), 
       labels = express(paste('beta[',1:ncol(sb), ']', sep='')),
       srt = 0,xpd = TRUE, cex=1, font=3) 
  abline(h=0.05, col='black',lwd=2)
  abline(h=0, col='black')
  text(x=ncol(sb)+3,y=0.08,labels='5%',adj=0, col='black',font=3)
  dev.off()
  
  # combine the results
  caseAll <- c('Cox-LASSO','GenF-f-pval','LogN-f-pval',
               'GenF-b-pval','LogN-b-pval','GenF-f-AIC',
               'LogN-f-AIC','GenF-b-AIC','LogN-b-AIC','GenF-f-BIC',
               'LogN-f-BIC','GenF-b-BIC','LogN-b-BIC','AFT-Bayes-LASSO')
  caseId <- c(13,7,1,10,4,8,2,11,5,9,3,12,6,14)
  len <- length(caseAll)
  
  for(mycase in c(1,2,4,6,8,10,12,14)){     #1:len
    load(paste('./genF/s1/c',caseId[mycase], sep=''))
    bw.dat$grp <- paste(bw.dat$grp, caseAll[mycase])
    if(mycase == 1) mat1 <- bw.dat
    if(mycase > 1) mat1 <- rbind(mat1, bw.dat)
  }
  bw.dat <- mat1
  
  simu_design <- 1
  
  mymain <- paste('under scenario 1 (30% censoring)')
  mains <- c('False Negatives','False Positives')
  k0 <- 1
  yadj <- c(-0.25, -1)
  
  for(k in c('N','P')){
    postscript(file=paste0('./box4rev_s',simu_design,'_',k,'.eps'),
               pointsize=10,width=7,height=3,horizontal=F,paper='special')  
    tmp <- bw.dat[which(substr(bw.dat$grp,7,7)==k),]
    tmp$grp <- factor(tmp$grp, levels=unique(tmp$grp))
    levels(tmp$grp) <- substr(levels(tmp$grp), 17, 100L)
    len <- length(levs <- levels(tmp$grp))
    par(las=1, mar=c(4.4,1.8,1.4,0)+.4,mgp=c(1.2,.3,0), tck=-0.02, 
        cex.axis=1, cex.lab=1, cex.main=1, xaxt='n') 
    boxplot(y~grp, data=tmp, ylab='Frequency', main=paste(mains[k0],mymain))
    text(c(1:len)+.3,rep(-0.25,len), labels = levs ,srt = 25,adj = c(1.1,1.1),xpd = TRUE, cex=1)
    k0 <- k0+1
    dev.off()
  }
  
  load('./genF/mysBeta2_GenF-b-pval-s2')
  sb3 <- mysBeta_gf
  
  falsePostive3 <- rowSums(sb3[,5:ncol(sb3)])
  falseNegatives3 <- 4-rowSums(sb3[,1:4])
  bw3.dat <- as.data.frame(c(falsePostive3, falseNegatives3)); names(bw3.dat) <- 'y'
  bw3.dat$grp <- as.factor(c(rep('False Positives',length(falsePostive3)),
                             rep('False Negatives',length(falseNegatives3))))
  mean(bw3.dat$y[bw3.dat$grp == 'False Positives']/16) 
  mean(bw3.dat$y[bw3.dat$grp == 'False Negatives']/4)

  ## manually
  bw_s2_GenF_b_pval <- bw3.dat
  levels(bw_s2_GenF_b_pval$grp) <- paste(
    levels(bw_s2_GenF_b_pval$grp), 
    "GenF-b-pval"
  )
  load('box_dat')
  
  for(simu_design in c(1,2)){
    mains <- c('False Negatives','False Positives')
    
    # for simulation 1
    if(simu_design == 1){
      bw.dat <- rbind(bw2.dat, bw_s1_genf_f_p,bw_s1_lnorm_f_p,bw_s1_genf_b_p, bw_s1_lnorm_b_p, bw_s1_genf_f,bw_s1_LogN_f_AIC,bw_s1_genf_b,bw_s1_LogN_b_AIC,bw_s1_GenF_f_BIC,bw_s1_LogN_f_BIC, bw_s1_genf_b_bic,bw_s1_LogN_b_BIC, bw_s1)
      mymain <- paste('for Simulation design 1 (30% censoring)')
    }
    
    # for simulation 2
    if(simu_design == 2){ 
      bw.dat <- rbind(bw_cox_lasso_s2.dat, bw_s2_GenF_f_pval, bw_s2_LogN_f_pval, bw_s2_GenF_b_pval, bw_s2_LogN_b_pval, bw_s2_GenF_f_AIC, bw_s2_LogN_f_AIC, bw_s2_GenF_b_AIC, bw_s2_LogN_b_AIC, bw_s2_GenF_f_BIC, bw_s2_LogN_f_BIC, bw_s2_GenF_b_BIC, bw_s2_LogN_b_BIC, bw_s2)
      mymain <- paste('for Simulation design 2 (70% censoring)')
    }
    
    load('mylevs0')
    levels(bw.dat$grp) <- mylevs0
    len <- length(levs <- levels(bw.dat$grp))
    
    mains <- c('False Negatives','False Positives')
    k0 <- 1
    yadj <- c(-0.25, -1)
    
    for(k in c('N','P')){
      postscript(file=paste0('./box4rev_s',simu_design,'_',k,'.eps'),
                 pointsize=10,width=7,height=3,horizontal=F,paper='special')  
      tmp <- bw.dat[which(substr(bw.dat$grp,7,7)==k),]
      tmp$grp <- factor(tmp$grp, levels=unique(tmp$grp))
      levels(tmp$grp) <- substr(levels(tmp$grp), 17, 100L)
      len <- length(levs <- levels(tmp$grp))
      par(las=1, mar=c(4.4,1.8,1.4,0)+.4,mgp=c(1.2,.3,0), 
          tck=-0.02, cex.axis=1, cex.lab=1, cex.main=1, xaxt='n') 
      boxplot(y~grp, data=tmp, ylab='Frequency', main=paste(mains[k0],mymain))
      text(c(1:len)+.3,rep(yadj[k0],len), labels = levs ,srt = 25,adj = c(1.1,1.1),xpd = TRUE, cex=1)
      k0 <- k0+1
      dev.off()
    }
    
  }
  
  levels(bw_s2$grp) <- paste(levels(bw_s2$grp), "AFT")
  levels(bw3.dat$grp) <- paste(levels(bw3.dat$grp), "Cox-LASSO")
  
  bw.dat <- rbind(bw_s2, bw3.dat)
  par(las=1, mar=c(1.4,2.2,.4,0)+.4,mgp=c(1.8,.3,0), tck=-0.01, cex.axis=1, cex.lab=1, cex.main=1) 
  boxplot(y~grp, data=bw.dat, ylab='Frequency')
  
}

if(runChunk == 3){
  ##----------- chunk 3: run in HPC for different repl
  # run by hierarchy
  library(R.matlab)
  library(flexsurv)
  library(survival)
  
  options(warn=-1) #disable warnings
  # argument repl indicates the bin number
  repl <- as.numeric(commandArgs(TRUE))  # 1-125
  
  load('realdat4genf')
  p <- ncol(tmp)-2
  load('out0')
  len <- length(indc <- 1:p)
  if(length(inds)>0) len <- length(indc <- indc[-inds])
  aic0 <- AIC(myfit)
  aics <- aic0-1
  
  if(all(indc!=repl)) stop('no need to proceed!')
  
  ptm <- proc.time()[3]
  q <- length(inds0 <- c(inds, repl))
  cat('q = ',q,':  ',sep='')
  cat('X',repl,' ',sep='')
  #if(!i%%10) cat('\n')
  
  myform0 <-  as.formula(paste("Surv(time, status) ~ ", paste(names(tmp)[inds0], collapse= "+")))
  fit0 <- flexsurvreg(myform0, data=tmp, dist='lnorm',method = "Nelder-Mead")
  coefs <- as.numeric(coef(fit0))
  init0 <- c(coefs[1],exp(coefs[2]),0, 1e-6, coefs[-c(1:2)])
  myfit <- update(fit0, inits=init0, data=tmp, dist='genf',method = "Nelder-Mead")
  aics <- AIC(myfit)
  
  cat('\n')
  cputime <- as.numeric(proc.time()[3]-ptm)
  cputime <- cputime/60
  
  cat('search steps done with CPUtime', round(cputime,3), 'minutes: completed!','\n')
  save(file=paste('out',q,'-',repl,sep=''),'myfit','cputime')
  
  
}


