NP_nst_transfo <- function(theta){

  return(c(log(theta[1]),                                              # sigmaF3
           log(theta[2]),                                              # sigmaF4
           log((1+theta[3])/(1-theta[3])),                             # rho
           log(theta[4]),                                              # sigmaR
           log(theta[5]),                                              # sigmaN
           log(theta[6]),                                              # sigmaP
           log(theta[7]),                                              # sigmaC
           log(theta[8]),                                              # q3
           log(theta[9]),                                              # q4
           log(theta[10]),                                             # q5
           log(theta[11]),                                             # q6
           log(theta[12]),                                             # q7
           log(theta[13]),                                             # q8
           log(theta[14])                                              # sigmaI
  ))
}
# END NP_nst_transfo
