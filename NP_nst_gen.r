NP_nst_gen <- function(theta,logF1=NULL,logN1=NULL,seedval=NULL){
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  require(MASS)
  sigmaF3   <- theta[1]
  sigmaF4   <- theta[2]
  rho       <- theta[3]
  sigmaR    <- theta[4]
  sigmaN    <- theta[5]
  sigmaP    <- theta[6]
  sigmaC    <- theta[7]
  q3        <- theta[8]
  q4        <- theta[9]
  q5        <- theta[10]
  q6        <- theta[11]
  q7        <- theta[12]
  q8        <- theta[13]
  sigmaI    <- theta[14] # p=14 for NP_nst
  
  age.F <- c(3:8,'9+')
  year.F <- 1967:2016
  AF <- length(age.F)
  TF <- length(year.F)
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
  
  if (!is.null(seedval)){set.seed(seedval)}
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate F ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### setup
  Sigmaxi <- matrix(NA_real_,AF,AF)
  Sigmaxi[1,1] <- sigmaF3^2
  Sigmaxi[1,-1] <- rho^(1:(AF-1))*sigmaF3*sigmaF4 # 1st row - [1,1]
  Sigmaxi[-1,1] <- rho^(1:(AF-1))*sigmaF4*sigmaF3 # 1st col - [1,1]
  Sigmaxi[-1,-1] <- rho^abs(outer(1:(AF-1),1:(AF-1),"-"))*sigmaF4^2
  
  logFat <- matrix(NA_real_,AF,TF)
  
  ### ini dist logF
  if (is.null(logF1)){ # needs fixed ini cond
    stop('Please supply logF1.')
  } else {
    logFat[,1] <- as.numeric(logF1) # fixed ini cond
  }
  
  ### dynamics for logF
  for (t in 2:TF){
    xit <- mvrnorm(n=1,mu=rep(0,AF),Sigma=Sigmaxi)
    logFat[,t] <- logFat[,t-1]+xit # RW
  }
  
  ### Fat and Zat
  Fat <- exp(logFat) # (AF x TF)
  Zat <- Mat+rbind(Fat,Fat[AF,]) # same dim as Mat (AN x TN)
  logZat <- log(Zat) # same dim as Mat (AN x TN)
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate N ####
  #/////////////////////////////////////////////////////////////////////////////
  
  ### setup
  logNat <- matrix(NA_real_,AN,TN)
  
  ### ini dist logN
  if (is.null(logN1)){ # needs fixed ini cond
    stop('Please supply logN1.')
  } else {
    logNat[,1] <- as.numeric(logN1) # fixed ini cond
  }
  
  Nat <- exp(logNat) # (AN x TN)
  
  ### dynamics for logN
  for (t in 2:TN){
    # RW for recruits
    logNat[1,t] <- logNat[1,t-1] + rnorm(1,0,sigmaR)
    Nat[1,t] <- exp(logNat[1,t])
    # survival of middle ages
    for (j in 2:(AN-1)){
      logNat[j,t] <- logNat[j-1,t-1]-Fat[j-1,t-1]-Mat[j-1,1] + rnorm(1,0,sigmaN)
      Nat[j,t] <- exp(logNat[j,t])
    }
    # survival of plus-group, Fat fixed at a=AF
    logNat[AN,t] <- log(Nat[AN-1,t-1]*exp(-Fat[AF,t-1]-Mat[AN-1,1]) +
                          + Nat[AN,t-1]*exp(-Fat[AF,t-1]-Mat[AN,1])) +
      + rnorm(1,0,sigmaP)
    Nat[AN,t] <- exp(logNat[AN,t])
  }
  

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
              'Mat'=Mat,'Zat'=Zat,
              'theta'=theta,'daysprop'=daysprop,
              'logF1'=logF1,'logN1'=logN1
              ))
}
# END NP_nst_gen
