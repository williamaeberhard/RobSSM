NP_st_weights <- function(obj,logF1=NULL,logN1=NULL){

  # tc[1] = meanlogF, phiF, sigmaF3, sigmaF4, rho
  # tc[2] = phiR, sigmaR             
  # tc[3] = meanlogN3, phiN, sigmaN
  # tc[4] = phiP, sigmaP          
  # tc[5] = sigmaC                   
  # tc[6] = q, sigmaI
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup ####
  #/////////////////////////////////////////////////////////////////////////////
  
  if (!is.null(logF1) | !is.null(logN1)){
    warning('logF1 and logN1 ignored, only st dist as ini implemented for now.')
  }
  
  require(mvtnorm) # for dmvnorm
  
  ### extract data and param from obj
  logCat <- obj$env$data$log_Cat
  logIat <- obj$env$data$log_Iat
  Mat <- obj$env$data$Mat
  daysprop <- obj$env$data$daysprop
  boundARcoef <- obj$env$data$boundARcoef
  tc <- obj$env$data$tc
  robcode <- obj$env$data$robcode
  ntrunc <- obj$env$data$ntrunc
  
  parfull <- obj$env$last.par.best
  whichrandeff <- obj$env$random # index in parfull corresponding to randeff
  
  theta.t <- parfull[-whichrandeff]
  theta <- NP_st_untransfo(theta.t)
  
  AF <- 7  # nrow(logFat) # age = 3, ..., 9+
  TF <- 50 # ncol(logFat) # t = 1967, ..., 2016
  AN <- 8  # nrow(logNat) # age = 3, ..., 10+
  TN <- 50 # ncol(logNat) # t = 1967, ..., 2016
  
  logFat <- matrix(parfull[names(parfull)=='log_Fat'],AF,TF)
  logNat <- matrix(parfull[names(parfull)=='log_Nat'],AN,TN)
  
  ### create rhoprime function according to robcode, loglog or SSH
  if (robcode==1){
    rhoprime <- function(x,tc){ # 1st deriv of loglog rho = weight
      # exp(x+tc)/(1+exp(x+tc))
      1-1/(1+exp(x+tc))
    }
  } else if (robcode==2){
    rhoprime <- function(x,tc){ # 1st deriv of SSH rho = weight
      ifelse(x>=(-tc),1,1/sqrt(1+((x+tc)/tc)^2))
    }
  } else {stop('robcode must be 1 (loglog) or 2 (SSH).')}
  
  ### extract design
  # meanlogF3 <- theta[1]
  # meanlogF4 <- theta[2]
  # meanlogF5 <- theta[3]
  # meanlogF6 <- theta[4]
  # meanlogF7 <- theta[5]
  # meanlogF8 <- theta[6]
  # meanlogF9 <- theta[7]
  phiF <-      theta[8]
  # sigmaF3 <-   theta[9]
  # sigmaF4 <-   theta[10]
  # rho <-       theta[11]
  # meanlogN3 <- theta[12]
  phiR <-      theta[13]
  sigmaR <-    theta[14]
  phiN <-      theta[15]
  sigmaN <-    theta[16]
  phiP <-      theta[17]
  sigmaP <-    theta[18]
  sigmaC <-    theta[19]
  # q3 <-        theta[20]
  # q4 <-        theta[21]
  # q5 <-        theta[22]
  # q6 <-        theta[23]
  # q7 <-        theta[24]
  # q8 <-        theta[25]
  sigmaI <-    theta[26]  # p=26
  
  AC <- nrow(logCat) # 8   # age = 3, ..., 10+
  TC <- ncol(logCat) # 49  # t = 1967, ..., 2015
  AI <- nrow(logIat) # 6   # age = 3, ..., 8+
  TI <- ncol(logIat) # 25  # t = 1992, ..., 2016
  t1992 <- TN-TI # time offset for variables ranging 1967-2016
  
  w1 <- double(TF) # multnorm of logFat
  w2 <- matrix(NA_real_,AN,TN) # dnorm of logNat
  w3 <- matrix(NA_real_,AC,TC) # dnorm of logCat
  w4 <- matrix(NA_real_,AI,TI) # dnorm of logIat
  
  Fat <- exp(logFat) #  same dim as logFat (AF x TF)
  Nat <- exp(logNat) # same dim as logNat (AN x TN)
  Zat <- Mat+rbind(Fat,Fat[AF,]) # same dim as Mat (AN x TN)
  logZat <- log(Zat) # same dim as Mat (AN x TN)
  
  ### compute st dist
  stdist <- NP_st_stdist(theta) # st dit as ini for both logF and logN
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### w1: weights on log F ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### ini
  # w1[1] <- 1 # starts at t=2, first weights fixed to 1
  w1[1] <- rhoprime(dmvnorm(x=logFat[,1],
                            mean=stdist$meanlogF,
                            sigma=stdist$varlogF,log=T),tc[1])
  # ^ st dist
  
  ### dynamics
  for (t in 2:TF){
    w1[t] <- rhoprime(dmvnorm(x=logFat[,t],
                              mean=(1-phiF)*stdist$meanlogF+phiF*logFat[,t-1],
                              sigma=stdist$Sigmaxi,log=T),tc[1])
    # ^ modified AR(1), st mean = meanlogF
  }
  
  w1.rescaled <- w1/max(w1) # rescale, relative to max weight
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### w2: weights on log N ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### ini
  w2[1,1] <- rhoprime(dnorm(x=logNat[1,1],
                            mean=stdist$meanlogN[1],
                            sd=sqrt(stdist$varlogN[1]),log=T),tc[2])
  for (a in 2:(AN-1)){
    w2[a,1] <- rhoprime(dnorm(x=logNat[a,1],
                              mean=stdist$meanlogN[a],
                              sd=sqrt(stdist$varlogN[a]),log=T),tc[3])
  }
  w2[AN,1] <- rhoprime(dnorm(x=logNat[AN,1],
                             mean=stdist$meanlogN[AN],
                             sd=sqrt(stdist$varlogN[AN]),log=T),tc[4])
  # ^ different tc for recruits, survival middle ages, survival plus-group
  
  ### dynamics
  for (t in 2:TF){
    w2[1,t] <- rhoprime(dnorm(x=logNat[1,t],
                              mean=(1-phiR)*stdist$meanlogN[1]+phiR*logNat[1,t-1],
                              sd=sigmaR,log=T),tc[2])
    for (a in 2:(AN-1)){
      mu_logNat <- phiN*(logNat[a-1,t-1]-Fat[a-1,t-1]-Mat[a-1,t-1])
      w2[a,t] <- rhoprime(dnorm(x=logNat[a,t],
                                mean=mu_logNat,
                                sd=sigmaN,log=T),tc[3])
    }
    mu_logNAt <- phiN*(logNat[AN-1,t-1]-Fat[AF,t-1]-Mat[AN-1,t-1]) +
      + phiP*(logNat[AN,t-1]-Fat[AF,t-1]-Mat[AN,t-1])
    w2[AN,t] <- rhoprime(dnorm(x=logNat[AN,t],
                               mean=mu_logNAt,
                               sd=sigmaP,log=T),tc[4])
  }
  
  w2.rescaled <- w2/max(w2) # rescale, relative to max weight
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### w3: weights on log C ####
  #/////////////////////////////////////////////////////////////////////////////
  
  for (t in 1:TC){
    for (a in 1:AF){ # AF=AC-1
      mu_logCat <- logFat[a,t]-logZat[a,t]+log(1-exp(-Zat[a,t]))+logNat[a,t]
      w3[a,t] <- rhoprime(dnorm(x=logCat[a,t],
                                mean=mu_logCat,
                                sd=sigmaC,log=T),tc[5])
    }
    mu_logCAt <- logFat[AF,t]-logZat[AC,t]+log(1-exp(-Zat[AC,t]))+logNat[AC,t]
    w3[AC,t] <- rhoprime(dnorm(x=logCat[AC,t],
                               mean=mu_logCAt,
                               sd=sigmaC,log=T),tc[5])
  }
  
  w3.rescaled <- w3/max(w3) # rescale, relative to max weight
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### w4: weights on log I ####
  #/////////////////////////////////////////////////////////////////////////////
  
  logq <- theta.t[20:25] # log(c(q3,q4,q5,q6,q7,q8))
  for (t in 1:TI){
    for (a in 1:AI){
      mu_logIat <- logq[a]-Zat[a,t+t1992]*daysprop+logNat[a,t+t1992]
      w4[a,t] <- rhoprime(dnorm(x=logIat[a,t],
                                mean=mu_logIat,
                                sd=sigmaI,log=T),tc[6])
    }
  }
  
  w4.rescaled <- w4/max(w4) # rescale, relative to max weight
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('w.logFat'=w1.rescaled,'w.logNat'=w2.rescaled,
              'w.logCat'=w3.rescaled,'w.logIat'=w4.rescaled,
              'w.logFat.unscaled'=w1,'w.logNat.unscaled'=w2,
              'w.logCat.unscaled'=w3,'w.logIat.unscaled'=w4))
}
# END NP_st_weights
