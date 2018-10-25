NP_nst_untransfo <- function(theta_t){
  
  return(c(exp(theta_t[1]),                                        # log_sigmaF3
           exp(theta_t[2]),                                        # log_sigmaF4
           (1-exp(-theta_t[3]))/(1+exp(-theta_t[3])),              # t_rho
           exp(theta_t[4]),                                        # log_sigmaR
           exp(theta_t[5]),                                        # log_sigmaN
           exp(theta_t[6]),                                        # log_sigmaP
           exp(theta_t[7]),                                        # log_sigmaC
           exp(theta_t[8]),                                        # log_q3
           exp(theta_t[9]),                                        # log_q4
           exp(theta_t[10]),                                       # log_q5
           exp(theta_t[11]),                                       # log_q6
           exp(theta_t[12]),                                       # log_q7
           exp(theta_t[13]),                                       # log_q8
           exp(theta_t[14])                                        # log_sigmaI
  ))
}
# END NP_nst_untransfo
