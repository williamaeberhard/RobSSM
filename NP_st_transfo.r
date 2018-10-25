NP_st_transfo <- function(theta,boundARcoef=0.95){

  return(c(theta[1],                                                 # meanlogF3
           theta[2],                                                 # meanlogF4
           theta[3],                                                 # meanlogF5
           theta[4],                                                 # meanlogF6
           theta[5],                                                 # meanlogF7
           theta[6],                                                 # meanlogF8
           theta[7],                                                 # meanlogF9
           log((boundARcoef+theta[8])/(boundARcoef-theta[8])),       # phiF
           log(theta[9]),                                            # sigmaF3
           log(theta[10]),                                           # sigmaF4
           log((1+theta[11])/(1-theta[11])),                         # rho
           theta[12],                                                # meanlogN3
           log((boundARcoef+theta[13])/(boundARcoef-theta[13])),     # phiR
           log(theta[14]),                                           # sigmaR
           log((1+theta[15])/(1-theta[15])),                         # phiN
           log(theta[16]),                                           # sigmaN
           log((boundARcoef+theta[17])/(boundARcoef-theta[17])),     # phiP
           log(theta[18]),                                           # sigmaP
           log(theta[19]),                                           # sigmaC
           log(theta[20]),                                           # q3
           log(theta[21]),                                           # q4
           log(theta[22]),                                           # q5
           log(theta[23]),                                           # q6
           log(theta[24]),                                           # q7
           log(theta[25]),                                           # q8
           log(theta[26])                                            # sigmaI
  ))
}
# END NP_st_transfo
