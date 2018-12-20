#/////////////////////////////////////////////////////////
#### RobSSM: Robust Estimation for State Space Models ####
#/////////////////////////////////////////////////////////

# RobSSM version 0.1, updated 2018-12-20
# Written by William H. Aeberhard <william.aeberhard@gmail.com>

# Disclaimer:
#  Users of these routines are cautioned that, while due care has been taken and
#  they are believed accurate, they have not been rigorously tested and their
#  use and results are solely the responsibilities of the user.

# Reference paper for methodology details:
#   Aeberhard, W. H., Cantoni, E., Field, C. Künsch, H. R., Mills Flemming, J.,
#   and Xu, X. (2018) Robust Estimation for General State Space Models. Submitted.
# Reference paper for SSM used for illustrationt:
#   Nielsen, A. and Berg, C. W. (2014). Estimation of time-varying selectivity
#   in stock assessments using state-space models. Fisheries Research 158, 96-101. 

# Two versions of SSM based on North Sea pollock stock assessment model:
#  * nst: non-stationary process, very close to original "Model D" (aka SAM) in
#         Nielsen and Berg (2014).
#  * st: stationary process, same observation equations as nst



#/////////////////////////////////////////////////////////////////////////////
#### nst: Simulate data, contaminate it, fit SSM by ML and robust methods ####
#/////////////////////////////////////////////////////////////////////////////

rm(list=ls())

### load libraries and functions
library(TMB)
compile('NP_nst.cpp') # compile TMB C++ template, to do only once
dyn.load(dynlib("NP_nst")) # load the compiled dynamic library

source('NP_nst_transfo.r') # transform theta for optimization over reals
source('NP_nst_untransfo.r') # re-transform back theta_t to original scale
source('NP_nst_gen.r') # generate sample according to nst SSM
source('NP_nst_nrcorrect.r') # Newton-Raphson correction step Fisher consistency
source('NP_nst_weights.r') # compute robustness weights given theta and randeff


### create design following simulation study in paper
# set up fishing mortality randeff
F1 <- c(0.266, 0.384, 0.357, 0.354, 0.313, 0.282, 0.314)
# ^ initial conditions fixed at predictions from initial run of SAM
logF1 <- log(as.numeric(F1)) # log(F) for year=1967, age=3,...,9+
age.F <- c(3:8,'9+') # age vector for F_{a,t}
year.F <- 1967:2016 # year vector
AF <- length(age.F) # nb of age classes "A" for F_{a,t}
TF <- length(year.F) # nb of time points "T"

# set up abundance randeff
N1 <- c(141236, 81733, 57375, 7179, 4931, 1159, 755, 690)
# ^ initial conditions fixed at predictions from initial run of SAM
logN1 <- log(N1) # log(N) for year=1967, age=3,...,10+
age.N <- c(3:9,'10+') # age vector for F_{a,t}
year.N <- 1967:2016 # year vector
AN <- length(age.N) # nb of age classes "A" for N_{a,t}
TN <- length(year.N) # same "T"

# set up dimensions and indices for catch response variable
age.C <- c(3:9,'10+') # age vector for C_{a,t}
year.C <- 1967:2015 # same missingness pattern as in North Sea pollock data
AC <- length(age.C) # number of age classes for C_{a,t}
TC <- length(year.C) # number of time points for C_{a,t}

# set up dimensions and indices for survey indices response variable
age.I <- c(3:7,'8+') # age vector for I_{a,t}
year.I <- 1992:2016 # same missingness pattern as in North Sea pollock data
AI <- length(age.I) # number of age classes for I_{a,t}
TI <- length(year.I) # number of time points for I_{a,t}
t1992 <- TN-TI # time offset for variables ranging 1967-2016

AC*TC+AI*TI # 542 = total nb of observations accounting for missingness

# set up fixed covariates
daysprop <- 0.6218853 # days for survey index, from North Sea pollock data
Mat <- matrix(0.2,AN,TN) # natural mortality, constant as in pollock data

# theta from North Sea pollock nst robust estimates
sigmaF3   <- 0.20   
sigmaF4   <- 0.16   
rho       <- 0.88   
sigmaR    <- 0.45   
sigmaN    <- 0.17   
sigmaP    <- 0.22   
sigmaC    <- 0.22   
q3        <- 5.1e-05
q4        <- 8.4e-05
q5        <- 6.3e-05
q6        <- 4.2e-05
q7        <- 2.8e-05
q8        <- 3.0e-05
sigmaI    <- 0.63   

theta <- c(sigmaF3,sigmaF4,rho,sigmaR,sigmaN,sigmaP,sigmaC,
           q3,q4,q5,q6,q7,q8,sigmaI)
names.theta <- c('sigmaF3','sigmaF4','rho','sigmaR','sigmaN','sigmaP','sigmaC',
                 paste0('q',3:8),'sigmaI')
names(theta) <- names.theta
p <- length(theta) # 14 for nst

# initial values for theta (transformed scale) and randeff (log scale)
theta.t.ini <- list()
theta.t.ini$log_sigmaF3 <-  0
theta.t.ini$log_sigmaF4 <-  0
theta.t.ini$t_rho       <-  1
theta.t.ini$log_sigmaR  <-  0
theta.t.ini$log_sigmaN  <-  0
theta.t.ini$log_sigmaP  <-  0
theta.t.ini$log_sigmaC  <-  0
theta.t.ini$log_q3      <- -5
theta.t.ini$log_q4      <- -5
theta.t.ini$log_q5      <- -5
theta.t.ini$log_q6      <- -5
theta.t.ini$log_q7      <- -5
theta.t.ini$log_q8      <- -5
theta.t.ini$log_sigmaI  <-  0

names.theta.t <- names(theta.t.ini)
theta.ini <- NP_nst_untransfo(unlist(theta.t.ini))
names(theta.ini) <- names.theta
theta.ini # intial values for theta on original scale

logFat.ini <- matrix(-1,AF,TF) # F_{a,t}
logNat.ini <- matrix(c(10,10,10,8,8,7,6,6),AN,TN) # N_{a,t}
# ^ initial values for all randeff, somewhat in ballpark

parlist.ini <- c(theta.t.ini,list('log_Fat'=logFat.ini,
                                  'log_Nat'=logNat.ini))
# ^ parameter list fed to TMB's MakeADFun

# indices for extraction of estimates and predicted randeff from TMB's sdreport
which.logFat <- (p+1):(p+AF*TF)
which.logNat <- (p+AF*TF+1):(p+AF*TF+AN*TN)
which.Fat <- (2*p+AF*TF+AN*TN+1):(2*p+2*AF*TF+AN*TN)
which.Nat <- (2*p+2*AF*TF+AN*TN+1):(2*p+2*AF*TF+2*AN*TN)

# settings for robust estimator
tc.ssh1 <- c(2.0, 3.0, 2.0, 2.0, 2.0, 4.0) # for 1st step
tc.ssh <- c(1.2, 2.3, 1.0, 1.3, 1.3, 2.9) # target tc
robcode <- 2 # 0=ML, 1=loglog, 2=SSH
H <- 1000 # nb MC samples in 1-step NR correction for robust est

# settings for contamination
whichcontamC <- 34 # t=34 in simulated Cat = year 2000
addC <- -4 # under-reported catches on log scale at single year


### simulate sample at the model
set.seed(1234) # reproducibility

dat <- NP_nst_gen(theta,logF1=logF1,logN1=logN1)
str(dat) # simulates both randeff and observations, returns everything in a list


### ML estimation, at model
datalist.ml.am <- list('log_Cat'=dat$logCat,'log_Iat'=dat$logIat,
                       'Mat'=Mat,'daysprop'=daysprop,
                       'tc'=rep(0,6),'robcode'=0)
# ^ data list fed to TMB's MakeADFun

obj.ml.am <- MakeADFun(data=datalist.ml.am,
                       parameters=parlist.ini,
                       random=c('log_Nat','log_Fat'),
                       DLL="NP_nst",silent=T)

# MLE: minimize (Laplace-approximated) marginal negloglik
system.time(opt.ml.am <- nlminb(start=obj.ml.am$par,obj=obj.ml.am$fn,gr=obj.ml.am$gr,
                                control=list(eval.max=1000,iter.max=1000)))
# ^ about 2 seconds on laptop with 2.3 GHz Intel Core i7
opt.ml.am
# ^ seems to have converged

# report theta estimates and randeff on original scale
system.time(rep.ml.am <- sdreport(obj.ml.am,bias.correct=F))
# ^ about 2 seconds on laptop with 2.3 GHz Intel Core i7

summ.rep <- summary(rep.ml.am)
theta.ml.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1]
theta.se.ml.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.ml.am <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.ml.am <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.ml.am <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.ml.am <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.ml.am <- matrix(summ.rep[which.Fat,1],AF,TF)
Fat.se.ml.am <- matrix(summ.rep[which.Fat,2],AF,TF)
Nat.ml.am <- matrix(summ.rep[which.Nat,1],AN,TN)
Nat.se.ml.am <- matrix(summ.rep[which.Nat,2],AN,TN)

cbind(theta,theta.ml.am,theta.se.ml.am)
# ^ true theta, MLE and se


### Robust estimation, at model

# Step 1/3: uncorrected robust estimator, minimize robustified marginal negloglik
# Step 2/3: eval corrected score based on initial uncorrected robust estimate
# Step 3/3: correct robust estimate by single NR step

# TODO:
# * robust estimation at model
# * compare ML and rob est and randeff (graph) at model
# * do the same unde contam
# * do all the above for st version


### contaminate simulated response
logCatcont <- dat$logCat
logCatcont[1,whichcontamC] <- pmax(logCatcont[1,whichcontamC]+addC,1)
# ^ contamination: underreport catch in year 2000, floor at 1 on log scale













### END RobSSM_Main
