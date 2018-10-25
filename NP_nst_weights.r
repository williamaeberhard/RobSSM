NP_nst_weights <- function(obj){

  # tc[1] = sigmaF3, sigmaF4, rho
  # tc[2] = sigmaR             
  # tc[3] = sigmaN
  # tc[4] = sigmaP          
  # tc[5] = sigmaC                   
  # tc[6] = q, sigmaI
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup ####
  #/////////////////////////////////////////////////////////////////////////////
  
  require(mvtnorm) # for dmvnorm
  
  ### extract data and param from obj
  logCat <- obj$env$data$log_Cat
  logIat <- obj$env$data$log_Iat
  Mat <- obj$env$data$Mat
  daysprop <- obj$env$data$daysprop
  tc <- obj$env$data$tc
  robcode <- obj$env$data$robcode

  parfull <- obj$env$last.par.best
  whichrandeff <- obj$env$random # index in parfull corresponding to randeff
  
  theta.t <- parfull[-whichrandeff]
  theta <- NP_nst_untransfo(theta.t)
  
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
  sigmaF3   <- theta[1]
  sigmaF4   <- theta[2]
  rho       <- theta[3]
  sigmaR    <- theta[4]
  sigmaN    <- theta[5]
  sigmaP    <- theta[6]
  sigmaC    <- theta[7]
  # q3        <- theta[8]
  # q4        <- theta[9]
  # q5        <- theta[10]
  # q6        <- theta[11]
  # q7        <- theta[12]
  # q8        <- theta[13]
  sigmaI    <- theta[14] # p=14 for NP_nst
  
  AC <- nrow(logCat) # 8   # age = 3, ..., 10+
  TC <- ncol(logCat) # 49  # t = 1967, ..., 2015
  AI <- nrow(logIat) # 6   # age = 3, ..., 8+
  TI <- ncol(logIat) # 25  # t = 1992, ..., 2016
  t1992 <- TN-TI # time offset for variables ranging 1967-2016
  
  w1 <- double(TF) # multnorm of logFat, starts at t=2
  w2 <- matrix(NA_real_,AN,TN) # dnorm of logNat, starts at t=2
  w3 <- matrix(NA_real_,AC,TC) # dnorm of logCat, starts at t=1
  w4 <- matrix(NA_real_,AI,TI) # dnorm of logIat, starts at t=1
  
  Fat <- exp(logFat) #  same dim as logFat (AF x TF)
  Nat <- exp(logNat) # same dim as logNat (AN x TN)
  Zat <- Mat+rbind(Fat,Fat[AF,]) # same dim as Mat (AN x TN)
  logZat <- log(Zat) # same dim as Mat (AN x TN)

  #/////////////////////////////////////////////////////////////////////////////
  #### w1: weights on log F ####
  #/////////////////////////////////////////////////////////////////////////////
  
  Sigmaxi <- matrix(NA_real_,AF,AF)
  Sigmaxi[1,1] <- sigmaF3^2
  Sigmaxi[1,-1] <- rho^(1:(AF-1))*sigmaF3*sigmaF4 # 1st row - [1,1]
  Sigmaxi[-1,1] <- rho^(1:(AF-1))*sigmaF4*sigmaF3 # 1st col - [1,1]
  Sigmaxi[-1,-1] <- rho^abs(outer(1:(AF-1),1:(AF-1),"-"))*sigmaF4^2
  
  w1[1] <- 1 # starts at t=2, first weight fixed to 1
  
  for (t in 2:TF){ # TN=TF
    ### w1: proc F
    w1[t] <- rhoprime(dmvnorm(x=logFat[,t],
                              mean=logFat[,t-1], # RW
                              sigma=Sigmaxi,log=T),tc[1])
  }
  
  w1.rescaled <- w1
  w1.rescaled[-1] <- w1[-1]/max(w1[-1]) # rescale, relative to max weight
  
  #/////////////////////////////////////////////////////////////////////////////
  #### w2: weights on log N ####
  #/////////////////////////////////////////////////////////////////////////////
  
  w2[,1] <- rep(1,AN) # starts at t=2, first weights fixed to 1
  
  for (t in 2:TF){ # TN=TF
    w2[1,t] <- rhoprime(dnorm(x=logNat[1,t],
                              mean=logNat[1,t-1], # RW
                              sd=sigmaR,log=T),tc[2])
    for (a in 2:(AN-1)){
      mu_logNat <- logNat[a-1,t-1] - Fat[a-1,t-1] - Mat[a-1,t-1]
      w2[a,t] <- rhoprime(dnorm(x=logNat[a,t],
                                mean=mu_logNat,
                                sd=sigmaN,log=T),tc[3])
    }
    mu_logNAt  <-  log(Nat[AN-1,t-1]*exp(-Fat[AF,t-1]-Mat[AN-1,t-1]) +
                         Nat[AN,t-1]*exp(-Fat[AF,t-1]-Mat[AN,t-1]))
    w2[AN,t] <- rhoprime(dnorm(x=logNat[AN,t],
                               mean=mu_logNAt,
                               sd=sigmaP,log=T),tc[4])
  }
  
  w2.rescaled <- w2
  w2.rescaled[,-1] <- w2[,-1]/max(w2[,-1]) # rescale, relative to max weight
  
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
  
  logq <- theta.t[8:13] # log(c(q3,q4,q5,q6,q7,q8))
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
# END NP_nst_weights
