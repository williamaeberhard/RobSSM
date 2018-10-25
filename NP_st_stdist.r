NP_st_stdist <- function(theta,ntrunc=20){
  # AF, AN and Mat hardcoded below
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup #### 
  #/////////////////////////////////////////////////////////////////////////////
  
  # meanlogF3 <- theta[1]
  # meanlogF4 <- theta[2]
  # meanlogF5 <- theta[3]
  # meanlogF6 <- theta[4]
  # meanlogF7 <- theta[5]
  # meanlogF8 <- theta[6]
  # meanlogF9 <- theta[7]
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
  # sigmaC    <- theta[19]
  # q3        <- theta[20]
  # q4        <- theta[21]
  # q5        <- theta[22]
  # q6        <- theta[23]
  # q7        <- theta[24]
  # q8        <- theta[25]
  # sigmaI    <- theta[26] # p=26 for NP_st
  
  AF <- 7 # length(c(3:8,'9+'))
  AN <- 8 # length(c(3:9,'10+'))
  Mat <- rep(0.2,AN) # constant, taken from North Sea Pollock design
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate F ####
  #/////////////////////////////////////////////////////////////////////////////
  
  Sigmaxi <- matrix(NA_real_,AF,AF)
  Sigmaxi[1,1] <- sigmaF3^2
  Sigmaxi[1,-1] <- rho^(1:(AF-1))*sigmaF3*sigmaF4 # 1st row - [1,1]
  Sigmaxi[-1,1] <- rho^(1:(AF-1))*sigmaF4*sigmaF3 # 1st col - [1,1]
  Sigmaxi[-1,-1] <- rho^abs(outer(1:(AF-1),1:(AF-1),"-"))*sigmaF4^2
  
  meanlogF <- as.numeric(theta[1:AF])
  varlogF <- Sigmaxi/(1-phiF^2) # st varcov logF
  meanF <- exp(meanlogF+diag(Sigmaxi)/(2*(1-phiF^2))) # st mean of log-normal F
  # varF <- exp(outer(meanlogF,meanlogF,'+') + 
  #               + outer(diag(Sigmaxi),diag(Sigmaxi),'+')/(2*(1-phiF^2)))*
  #   (exp(Sigmaxi/(1-phiF^2))-1) # st varcov of log-normal F, not necessary
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Generate N ####
  #/////////////////////////////////////////////////////////////////////////////
    
  meanlogN <- double(AN)
  varlogN <- meanlogN # vector only
  
  meanlogN[1] <- meanlogN3 # by construction
  varlogN[1] <- sigmaR^2/(1-phiR^2) # st AR(1) for recruits
  
  for (a in 2:(AN-1)){ # ini, survival of middle ages
    meanlogN[a] <- phiN^(a-1)*meanlogN3 - sum(phiN^(1:(a-1))*
                                                (meanF[(a-1):1]+Mat[(a-1):1]))
    # ^ Mat constant through time
    varlogN[a] <- phiN^(2*a-2)*varlogN[1] +
      + sigmaN^2*(1-phiN^(2*(a-1)))/(1-phiN^2) +
      # ^ finite geometric sum = sigmaN^2*sum(phiN^(2*(0:(a-2))))
      + sum(phiN^outer(1:(a-1),1:(a-1),'+')*
              exp(outer(meanlogF[(a-1):1],meanlogF[(a-1):1],'+') + 
                    + outer(diag(Sigmaxi)[(a-1):1],
                            diag(Sigmaxi)[(a-1):1],'+')/(2*(1-phiF^2)))*
              (exp(phiF^abs(outer(1:(a-1),1:(a-1),"-"))*
                     Sigmaxi[(a-1):1,(a-1):1]/(1-phiF^2))-1))
    # ^ correct with Fat kept random!
  }
  
  meanlogN[AN] <- (phiN^(AN-1)*meanlogN3 + # ini plus group
                     -sum(phiN^(1:(AN-1))*(meanF[(AN-1):1]+Mat[(AN-1):1])) +
                     -phiP*(meanF[AF]+Mat[AN]))/(1-phiP) # plus group
  
  uV1 <- 0
  for (i in 1:ntrunc){
    for (j in 1:ntrunc){
      uV1 <- uV1 + phiP^(i+j)*(exp(phiF^abs(i-j)*varlogF[AF,AF])-1) # AN
    }
  }
  uV1 <- uV1*exp(2*meanlogF[AF]+varlogF[AF,AF]) # AN
  
  uV2 <- 0
  for (k in 1:(AN-1)){
    for (l in 1:(AN-1)){
      for (i in 0:ntrunc){
        for (j in 0:ntrunc){
          uV2 <- uV2 + phiP^(i+j)*phiN^(k+l)*
            exp(meanlogF[AN-k]+meanlogF[AN-l]+sum(diag(varlogF)[c(AN-k,AN-l)])/2)*
            (exp(phiF^abs(j+k-i-l)*varlogF[AN-k,AN-l])-1)
        }
      }
    }
  }
  
  uv3.alt <- 0
  for (j in 1:ntrunc){
    for (i in 0:ntrunc){
      for (k in 1:(AN-1)){
        uv3.alt <- uv3.alt + phiP^i*phiP^j*phiN^k*
          exp(meanlogF[AN-k]+meanlogF[AF]+(varlogF[AN-k,AN-k]+varlogF[AF,AF])/2)*
          (exp(phiF^abs(j-k-i)*varlogF[AN-k,AF])-1)
      }
    }
  }
  
  varlogN[AN] <- sigmaP^2/(1-phiP^2) + 
    + sigmaN^2*phiN^2*(1-phiN^(2*(AN-2)))/(1-phiN^2)/(1-phiP^2) +
    # ^ finite geometric sum = sigmaN^2*sum(phiN^(2*(1:(AN-2))))/(1-phiP^2)
    + sigmaR^2*phiN^(2*(AN-1))/(phiR-phiP)^2*
    (phiR^2/(1-phiR^2)-2*phiR*phiP/(1-phiR*phiP)+phiP^2/(1-phiP^2)) +
    + uV1 + uV2 + 2*uv3.alt # factor 2 is correct!
  
  meanN <- exp(meanlogN+varlogN/2) # log-normal mean, simplest approx
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  return(list('theta'=theta,'Sigmaxi'=Sigmaxi,
              'meanlogF'=meanlogF,'varlogF'=varlogF,'meanF'=meanF,
              'meanlogN'=meanlogN,'varlogN'=varlogN,'meanN'=meanN
  ))
}
# END NP_st_stdist
