NP_st_untransfo <- function(theta_t,boundARcoef=0.95){

  return(c(theta_t[1],                                             # meanlogF3
           theta_t[2],                                             # meanlogF4
           theta_t[3],                                             # meanlogF5
           theta_t[4],                                             # meanlogF6
           theta_t[5],                                             # meanlogF7
           theta_t[6],                                             # meanlogF8
           theta_t[7],                                             # meanlogF9
           boundARcoef*(1-exp(-theta_t[8]))/(1+exp(-theta_t[8])),  # t_phiF
           exp(theta_t[9]),                                        # log_sigmaF3
           exp(theta_t[10]),                                       # log_sigmaF4
           (1-exp(-theta_t[11]))/(1+exp(-theta_t[11])),            # t_rho
           theta_t[12],                                            # meanlogN3
           boundARcoef*(1-exp(-theta_t[13]))/(1+exp(-theta_t[13])),# t_phiR
           exp(theta_t[14]),                                       # log_sigmaR
           (1-exp(-theta_t[15]))/(1+exp(-theta_t[15])),            # t_phiN
           exp(theta_t[16]),                                       # log_sigmaN
           boundARcoef*(1-exp(-theta_t[17]))/(1+exp(-theta_t[17])),# t_phiP
           exp(theta_t[18]),                                       # log_sigmaP
           exp(theta_t[19]),                                       # log_sigmaC
           exp(theta_t[20]),                                       # log_q3
           exp(theta_t[21]),                                       # log_q4
           exp(theta_t[22]),                                       # log_q5
           exp(theta_t[23]),                                       # log_q6
           exp(theta_t[24]),                                       # log_q7
           exp(theta_t[25]),                                       # log_q8
           exp(theta_t[26])                                        # log_sigmaI
  ))
}
# END NP_st_untransfo
