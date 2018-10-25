NP_st_gen <- function(theta,ntrunc=20,logF1=NULL,logN1=NULL,seedval=NULL,
                      stdist=NULL){
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  require(MASS)
  meanlogF3 <- theta[1]
  meanlogF4 <- theta[2]
  meanlogF5 <- theta[3]
  meanlogF6 <- theta[4]
  meanlogF7 <- theta[5]
  meanlogF8 <- theta[6]
  meanlogF9 <- theta[7]
  phiF      <- theta[8]
  sigmaF3   <- theta[9]
  sigmaF4   <- theta[10]
  rho       <- theta[11]
  meanlogN3 <- theta[12]
  phiR      <- theta[13]
  sigmaR    <- theta[14]
  phiN      <- theta[15]
  sigmaN    <- theta[16]
  phiP      <- theta[17]
  sigmaP    <- theta[18]
  sigmaC    <- theta[19]
  q3        <- theta[20]
  q4        <- theta[21]
  q5        <- theta[22]
  q6        <- theta[23]
  q7        <- theta[24]
  q8        <- theta[25]
  sigmaI    <- theta[26] # p=26 for NP_st
  
  age.F <- c(3:8,'9+')
  year.F <- 1967:2016
  AF <- length(age.F) # 7
  TF <- length(year.F) # 50
  
  age.N <- c(3:9,'10+')
  year.N <- 1967:2016
  AN <- length(age.N) # 8
  TN <- length(year.N) # 50
  
  age.C <- c(3:9,'10+')
  year.C <- 1967:2015
  AC <- length(age.C) # number of age classes for catch # 8
  TC <- length(year.C) # number of time points for catch # 49
  
  age.I <- c(3:7,'8+')
  year.I <- 1992:2016
  AI <- length(age.I) # number of age classes for survey index # 6
  TI <- length(year.I) # number of time points for survey index # 25
  
  t1992 <- TN-TI # time offset for variables ranging 1967-2016
  daysprop <- 0.6218853 # from SaitheNS data
  Mat <- matrix(0.2,AN,TN) # from SaitheNS data
  
  if (is.null(stdist)){
    stdist <- NP_st_stdist(theta) # st dit as ini for both logF and logN
  }
  
  if (!is.null(seedval)){set.seed(seedval)}
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate F ####
  #/////////////////////////////////////////////////////////////////////////////
  
  logFat <- matrix(NA_real_,AF,TF)
  
  ### ini dist logF
  if (!is.null(logF1)){ # then fixed ini cond rather than st dist
    warning('Fixed ini cond for F, using supplied logF1.')
    logFat[,1] <- as.numeric(logF1) # fixed ini cond
  } else {
    logFat[,1] <- mvrnorm(n=1,mu=stdist$meanlogF,Sigma=stdist$varlogF) # st dist
  }
  
  ### dynamics for logF  
  for (t in 2:TF){
    xit <- mvrnorm(n=1,mu=rep(0,AF),Sigma=stdist$Sigmaxi)
    logFat[,t] <- (1-phiF)*stdist$meanlogF+phiF*logFat[,t-1]+xit # modified AR(1)
  }
  
  ### Fat and Zat
  Fat <- exp(logFat) # (AF x TF)
  Zat <- Mat+rbind(Fat,Fat[AF,]) # same dim as Mat (AN x TN)
  logZat <- log(Zat) # same dim as Mat (AN x TN)
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate N ####
  #/////////////////////////////////////////////////////////////////////////////
  
  logNat <- matrix(NA_real_,AN,TN)
  
  ### ini dist logN
  if (!is.null(logN1)){ # then fixed ini cond rather than st dist
    warning('Fixed ini cond for N, using supplied logN1.')
    logNat[,1] <- as.numeric(logN1) # fixed ini cond
  } else { # then use st dist
    # for (a in 1:AN){
    #   logNat[a,1] <- rnorm(n=1,mean=stdist$meanlogN[a],
    #                        sd=sqrt(stdist$varlogN[a]))
    # }
    logNat[,1] <- rnorm(n=AN,mean=stdist$meanlogN,sd=sqrt(stdist$varlogN))
    # ^ indep across ages
  }
  
  ### dynamics for logN
  for (t in 2:TN){
    # AR(1) for recruits
    logNat[1,t] <- (1-phiR)*stdist$meanlogN[1] + phiR*logNat[1,t-1] +
      + rnorm(1,0,sigmaR)
    # Nat[1,t] <- exp(logNat[1,t])
    # survival of middle ages
    for (j in 2:(AN-1)){
      logNat[j,t] <- phiN*(logNat[j-1,t-1]-Fat[j-1,t-1]-Mat[j-1,1])+
        + rnorm(1,0,sigmaN)
      # Nat[j,t] <- exp(logNat[j,t])
    }
    # survival of plus-group, Fat fixed at a=AF
    logNat[AN,t] <- phiN*(logNat[AN-1,t-1]-Fat[AF,t-1]-Mat[AN-1,1]) +
      + phiP*(logNat[AN,t-1]-Fat[AF,t-1]-Mat[AN,1]) + rnorm(1,0,sigmaP)
    # Nat[AN,t] <- exp(logNat[AN,t])
  }
  
  Nat <- exp(logNat) # (AN x TN)
  

  #/////////////////////////////////////////////////////////////////////////////
  #### Generate C ####
  #/////////////////////////////////////////////////////////////////////////////
  
  logCat <- matrix(NA_real_,AC,TC)
  for (t in 1:TC){
    for (a in 1:AF){ # AF=AC-1
      logCat[a,t] <- logFat[a,t]-logZat[a,t]+
        log(1-exp(-Zat[a,t]))+logNat[a,t] + rnorm(1,0,sigmaC)
    }
    logCat[AC,t] <- logFat[AF,t]-logZat[AC,t]+
      log(1-exp(-Zat[AC,t]))+logNat[AC,t] + rnorm(1,0,sigmaC)
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate I ####
  #/////////////////////////////////////////////////////////////////////////////
  
  logIat <- matrix(NA_real_,AI,TI)
  logq <- log(c(q3,q4,q5,q6,q7,q8))
  for (t in 1:TI){
    for (a in 1:AI){
      logIat[a,t] <- logq[a]-Zat[a,t+t1992]*daysprop+logNat[a,t+t1992] +
        + rnorm(1,0,sigmaI)
    }
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('Fat'=Fat,'logFat'=logFat,
              'Nat'=Nat,'logNat'=logNat,
              'logCat'=logCat,'logIat'=logIat,
              'Mat'=Mat,'Zat'=Zat,'daysprop'=daysprop,
              'theta'=theta,
              'meanlogF'=stdist$meanlogF,'varlogF'=stdist$varlogF,
              'meanF'=stdist$meanF,
              'meanlogN'=stdist$meanlogN,'varlogN'=stdist$varlogN,
              'meanN'=stdist$meanN
              ))
}
# END NP_st_gen
