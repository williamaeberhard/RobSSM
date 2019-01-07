#/////////////////////////////////////////////////////////
#### RobSSM: Robust Estimation for State Space Models ####
#/////////////////////////////////////////////////////////

# RobSSM version 0.1, updated 2019-01-07
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
set.seed(1234) # for reproducibility

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
Fat.ml.am <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.ml.am <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

cbind(theta,theta.ml.am,theta.se.ml.am)
# ^ true theta, MLE and se


### Robust estimation, at model

# Step 1/3: uncorrected robust estimator, minimize robustified marginal negloglik
datalist.rob.am <- list('log_Cat'=dat$logCat,'log_Iat'=dat$logIat,
                        'Mat'=Mat,'daysprop'=daysprop,
                        'tc'=tc.ssh,'robcode'=robcode)
parlist.robu.am <- c(opt.ml.am$par,list('log_Fat'=logFat.ini,
                                        'log_Nat'=logNat.ini))
obj.robu.am <- MakeADFun(data=datalist.rob.am,
                         parameters=parlist.robu.am,
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_nst",silent=T)
system.time(opt.robu.am <- nlminb(start=obj.robu.am$par,obj=obj.robu.am$fn,
                                  gr=obj.robu.am$gr,
                                  control=list(eval.max=1000,iter.max=1000)))
# ^ about 2s on laptop
opt.robu.am
# ^ seems to have converged

system.time(rep.robu.am <- sdreport(obj.robu.am,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robu.am) # overwrite summ.rep, big object
logFat.robu.am <- matrix(summ.rep[which.logFat,1],AF,TF)
logNat.robu.am <- matrix(summ.rep[which.logNat,1],AN,TN)
# ^ predicted randeff based on uncorrected robust estimate, only for initial
#   values below


# Step 2/3 and 3/3: evaluate corrected score based on initial uncorrected robust
#                   estimate and correct it by single NR step

system.time(opt.robc.am <- NP_nst_nrcorrect(theta.t=opt.robu.am$par,
                                            obj=obj.robu.am,H=H,maxit=1))
# ^ XXX seconds with H=1000 # TODO
str(opt.robc.am,1)
# ^ list with corrected robust gradient and (finite diff numerical) Hessian


# predict randeff based on corrected theta
parlist.robc <- c(opt.robc.am$theta.t,list('log_Fat'=logFat.robu.am,
                                           'log_Nat'=logNat.robu.am))
obj.robc.am <- MakeADFun(data=datalist.rob.am, # same data
                         parameters=parlist.robc, # updated par
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_nst",silent=T)
invisible(obj.robc.am$fn())
# ^ no outer optim, just predict randeff with inner optim (Laplace approx)

system.time(rep.robc.am <- sdreport(obj.robc.am,bias.correct=F))
# ^ about 2 s on laptop

summ.rep <- summary(rep.robc.am) # overwrite summ.rep, big object
theta.robc.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1] # opt.robc.am$theta
theta.se.robc.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.robc.am <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.robc.am <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.robc.am <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.robc.am <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.robc.am <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.robc.am <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

w.rob.am <- NP_nst_weights(obj.robc.am) # robustness weights
w.rob.am <- w.rob.am[5:8] # no need to rescale
str(w.rob.am)
# ^ list of 4 matrices, one for each loglik contribution: F_{a,t}, N_{a,t},
#   C_{a,t}, and I_{a,t}


# compare ML and robust estimates at the model
cbind(theta,theta.ml.am,theta.se.ml.am,theta.robc.am,theta.se.robc.am)
# ^ true theta, MLE and its se, robust and its se

par(mfrow=c(5,3))
for (j in 1:p){
  plot(1:2,c(theta.ml.am[j],theta.robc.am[j]),pch=19,
       xlim=c(0.7,2.3),ylim=theta[j]*c(0.5,1.5),xaxt='n',xlab='',ylab='',
       main=names.theta[j])
  arrows(x0=1:2,y0=c(theta.ml.am[j]-2*theta.se.ml.am[j],
                     theta.robc.am[j]-2*theta.se.robc.am[j]),
         x1=1:2,y1=c(theta.ml.am[j]+2*theta.se.ml.am[j],
                     theta.robc.am[j]+2*theta.se.robc.am[j]),
         angle=90,code=3,length=0.05)
  # ^ error bars as rough 95% CI, just to give a sense of variability
  axis(1,at=1:2,labels=c('ML','Robust'),lwd=0,lwd.ticks=1)
  abline(h=theta[j],col='red')
}
par(mfrow=c(1,1))
# ^ MLE and robust est generally nearly identical at the model


# compare ML and robust predicted randeff at the model

lb.Fat.ml.am <- exp(logFat.ml.am-2*logFat.se.ml.am)
ub.Fat.ml.am <- exp(logFat.ml.am+2*logFat.se.ml.am)
lb.Nat.ml.am <- exp(logNat.ml.am-2*logNat.se.ml.am)
ub.Nat.ml.am <- exp(logNat.ml.am+2*logNat.se.ml.am)
# ^ exponentiate bounds back on original scale

lb.Fat.robc.am <- exp(logFat.robc.am-2*logFat.se.robc.am)
ub.Fat.robc.am <- exp(logFat.robc.am+2*logFat.se.robc.am)
lb.Nat.robc.am <- exp(logNat.robc.am-2*logNat.se.robc.am)
ub.Nat.robc.am <- exp(logNat.robc.am+2*logNat.se.robc.am)
# ^ exponentiate bounds back on original scale

yearvec <- rep(NA,TF)
every5 <- (0:9)*5+4
yearvec[every5] <- year.F[every5] # label every 5 years

col.pollockages <- c('#da5e16','#7a13be','#f00000','#000cf0',
                     '#00850b','#c70074','#00c1d6','#6E7B8B')
col.pollockages.env <- paste0(col.pollockages,'30') # semi-transparent

names.age.F <- paste0('a = ',age.F)
names.age.N <- paste0('a = ',age.N)
names.age.I <- paste0('a = ',age.I)
names.age.C <- names.age.N



par(mfrow=c(2,2))
# nst Fat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} ML predictions, at model',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.ml.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.ml.am[j,],ub.Fat.ml.am[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend('topright',names.age.F,lty=1,bty='n',
       col=col.pollockages[1:AF],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# nst Fat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} robust predictions, at model',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.robc.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.robc.am[j,],ub.Fat.robc.am[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
# nst Nat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} ML predictions, at model',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.ml.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.ml.am[j,],ub.Nat.ml.am[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend(x='topright',names.age.N,lty=1,bty='n',
       col=col.pollockages[1:AN],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# nst Nat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} robust predictions, at model',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.robc.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.robc.am[j,],ub.Nat.robc.am[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
par(mfrow=c(1,1))


# investigate robustness weights at the model
lapply(w.rob.am,function(x){summary(as.numeric(x))})
# ^ nearly all weights are 1 as expected

unlist(lapply(w.rob.am,function(x){which(x<0.8,arr.ind=T)}))
# ^ only one really downweighted observation, given arbitrary threshold 0.8

dat$logCat[,11:15]
# ^ downweighted obs (6,13) may be atypically high given neighborhood


### contaminate simulated response
logCatcont <- dat$logCat
logCatcont[1,whichcontamC] <- pmax(logCatcont[1,whichcontamC]+addC,1)
# ^ contamination: underreport catch for recruits in year 2000, with floor at 1
#   on log scale to keep realistic values


### ML estimation, under contam
datalist.ml.uc <- list('log_Cat'=logCatcont,'log_Iat'=dat$logIat, # contaminated
                       'Mat'=Mat,'daysprop'=daysprop,
                       'tc'=rep(0,6),'robcode'=0)
# ^ data list fed to TMB's MakeADFun

obj.ml.uc <- MakeADFun(data=datalist.ml.uc, # contam data
                       parameters=parlist.ini, # same as at model
                       random=c('log_Nat','log_Fat'),
                       DLL="NP_nst",silent=T)

# MLE: minimize (Laplace-approximated) marginal negloglik
system.time(opt.ml.uc <- nlminb(start=obj.ml.uc$par,obj=obj.ml.uc$fn,gr=obj.ml.uc$gr,
                                control=list(eval.max=1000,iter.max=1000)))
# ^ about 2s
opt.ml.uc
# ^ seems to have converged

# report theta estimates and randeff on original scale
system.time(rep.ml.uc <- sdreport(obj.ml.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.ml.uc)
theta.ml.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1]
theta.se.ml.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.ml.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.ml.uc <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.ml.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.ml.uc <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.ml.uc <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.ml.uc <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

cbind(theta,theta.ml.am,theta.ml.uc)
# ^ true theta, MLE at model and MLE under contam: sigmaF3-4 most biased


### Robust estimation, under contam

# Step 1/3: uncorrected robust estimator, minimize robustified marginal negloglik
datalist.rob.uc <- list('log_Cat'=logCatcont,'log_Iat'=dat$logIat, # contam
                        'Mat'=Mat,'daysprop'=daysprop,
                        'tc'=tc.ssh,'robcode'=robcode)
parlist.robu.uc <- c(opt.ml.uc$par,list('log_Fat'=logFat.ini,
                                        'log_Nat'=logNat.ini))
obj.robu.uc <- MakeADFun(data=datalist.rob.uc,
                         parameters=parlist.robu.uc,
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_nst",silent=T)
system.time(opt.robu.uc <- nlminb(start=obj.robu.uc$par,obj=obj.robu.uc$fn,
                                  gr=obj.robu.uc$gr,
                                  control=list(eval.max=1000,iter.max=1000)))
# ^ about 3s
opt.robu.uc
# ^ seems to have converged

system.time(rep.robu.uc <- sdreport(obj.robu.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robu.uc) # overwrite summ.rep, big object
logFat.robu.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logNat.robu.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
# ^ predicted randeff based on uncorrected robust estimate, only for initial
#   values below


# Step 2/3 and 3/3: evaluate corrected score based on initial uncorrected robust
#                   estimate and correct it by single NR step

system.time(opt.robc.uc <- NP_nst_nrcorrect(theta.t=opt.robu.uc$par,
                                            obj=obj.robu.uc,H=H,maxit=1))
# ^ XXX seconds with H=1000 # TODO
str(opt.robc.uc,1)
# ^ list with corrected robust gradient and (finite diff numerical) Hessian


# predict randeff based on corrected theta
parlist.robc <- c(opt.robc.uc$theta.t,list('log_Fat'=logFat.robu.uc,
                                           'log_Nat'=logNat.robu.uc))
obj.robc.uc <- MakeADFun(data=datalist.rob.uc, # same data
                         parameters=parlist.robc, # updated par
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_nst",silent=T)
invisible(obj.robc.uc$fn())
# ^ no outer optim, just predict randeff with inner optim (Laplace approx)

system.time(rep.robc.uc <- sdreport(obj.robc.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robc.uc) # overwrite summ.rep, big object
theta.robc.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1] # opt.robc.uc$theta
theta.se.robc.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.robc.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.robc.uc <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.robc.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.robc.uc <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.robc.uc <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.robc.uc <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

w.rob.uc <- NP_nst_weights(obj.robc.uc) # robustness weights
w.rob.uc <- w.rob.uc[5:8] # no need to rescale
str(w.rob.uc)
# ^ list of 4 matrices, one for each loglik contribution: F_{a,t}, N_{a,t},
#   C_{a,t}, and I_{a,t}


# compare ML and robust estimates under contam
cbind(theta,theta.ml.uc,theta.se.ml.uc,theta.robc.uc,theta.se.robc.uc)
# ^ true theta, MLE and its se, robust and its se

cbind(abs(theta.ml.am-theta.ml.uc)/theta,abs(theta.robc.am-theta.robc.uc)/theta)
# ^ est diff between at model and under contam, relative to true theta: MLE
#   changed a lot for a few param, robust remains globally stable

par(mfrow=c(5,3))
for (j in 1:p){
  plot(1:2,c(theta.ml.uc[j],theta.robc.uc[j]),pch=19,
       xlim=c(0.7,2.3),ylim=theta[j]*c(0.5,1.5),xaxt='n',xlab='',ylab='',
       main=names.theta[j])
  arrows(x0=1:2,y0=c(theta.ml.uc[j]-2*theta.se.ml.uc[j],
                     theta.robc.uc[j]-2*theta.se.robc.uc[j]),
         x1=1:2,y1=c(theta.ml.uc[j]+2*theta.se.ml.uc[j],
                     theta.robc.uc[j]+2*theta.se.robc.uc[j]),
         angle=90,code=3,length=0.05)
  # ^ error bars as rough 95% CI, just to give a sense of variability
  axis(1,at=1:2,labels=c('ML','Robust'),lwd=0,lwd.ticks=1)
  abline(h=theta[j],col='red')
}
par(mfrow=c(1,1))
# ^ robust est remains close to true value, MLE biased now


# compare ML and robust predicted randeff at the model
lb.Fat.ml.uc <- exp(logFat.ml.uc-2*logFat.se.ml.uc)
ub.Fat.ml.uc <- exp(logFat.ml.uc+2*logFat.se.ml.uc)
lb.Nat.ml.uc <- exp(logNat.ml.uc-2*logNat.se.ml.uc)
ub.Nat.ml.uc <- exp(logNat.ml.uc+2*logNat.se.ml.uc)
# ^ exponentiate bounds back on original scale

lb.Fat.robc.uc <- exp(logFat.robc.uc-2*logFat.se.robc.uc)
ub.Fat.robc.uc <- exp(logFat.robc.uc+2*logFat.se.robc.uc)
lb.Nat.robc.uc <- exp(logNat.robc.uc-2*logNat.se.robc.uc)
ub.Nat.robc.uc <- exp(logNat.robc.uc+2*logNat.se.robc.uc)
# ^ exponentiate bounds back on original scale


par(mfrow=c(2,2))
# nst Fat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} ML predictions, under contam',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.ml.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.ml.uc[j,],ub.Fat.ml.uc[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend('topright',names.age.F,lty=1,bty='n',
       col=col.pollockages[1:AF],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# nst Fat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} robust predictions, under contam',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.robc.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.robc.uc[j,],ub.Fat.robc.uc[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
# nst Nat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} ML predictions, under contam',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.ml.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.ml.uc[j,],ub.Nat.ml.uc[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend(x='topright',names.age.N,lty=1,bty='n',
       col=col.pollockages[1:AN],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# nst Nat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} robust predictions, under contam',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.robc.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.robc.uc[j,],ub.Nat.robc.uc[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
par(mfrow=c(1,1))
# ^ ML predictions now biased around contminated year 2000 and trajectory of
#   recruits very jagged overall, robust predictions remain stable as expected


# investigate robustness weights under contamintion
lapply(w.rob.uc,function(x){summary(as.numeric(x))})
# ^ most weights are 1, but some very close to 0

lapply(w.rob.uc,function(x){which(x<0.8,arr.ind=T)})
# ^ one additional heavily downweighted commercial catch observation

logCatcont[,32:36]
# ^ downweighted obs (1,34) is exactly contaminated point in year 2000
dat$logCat[,32:36] # (1,34) was originally around 8



#/////////////////////////////////////////////////////////////////////////////
#### st: Simulate data, contaminate it, fit SSM by ML and robust methods ####
#/////////////////////////////////////////////////////////////////////////////

# here!!!

rm(list=ls())

### load libraries and functions
library(TMB)
compile('NP_st.cpp') # compile TMB C++ template, to do only once
dyn.load(dynlib("NP_st")) # load the compiled dynamic library

source('NP_st_transfo.r') # transform theta for optimization over reals
source('NP_st_untransfo.r') # re-transform back theta_t to original scale
source('NP_st_gen.r') # generate sample according to st SSM
source('NP_st_nrcorrect.r') # Newton-Raphson correction step Fisher consistency
source('NP_st_weights.r') # compute robustness weights given theta and randeff


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

# theta from North Sea pollock st robust estimates
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
p <- length(theta) # 14 for st

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
theta.ini <- NP_st_untransfo(unlist(theta.t.ini))
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
set.seed(1234) # for reproducibility

dat <- NP_st_gen(theta,logF1=logF1,logN1=logN1)
str(dat) # simulates both randeff and observations, returns everything in a list


### ML estimation, at model
datalist.ml.am <- list('log_Cat'=dat$logCat,'log_Iat'=dat$logIat,
                       'Mat'=Mat,'daysprop'=daysprop,
                       'tc'=rep(0,6),'robcode'=0)
# ^ data list fed to TMB's MakeADFun

obj.ml.am <- MakeADFun(data=datalist.ml.am,
                       parameters=parlist.ini,
                       random=c('log_Nat','log_Fat'),
                       DLL="NP_st",silent=T)

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
Fat.ml.am <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.ml.am <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

cbind(theta,theta.ml.am,theta.se.ml.am)
# ^ true theta, MLE and se


### Robust estimation, at model

# Step 1/3: uncorrected robust estimator, minimize robustified marginal negloglik
datalist.rob.am <- list('log_Cat'=dat$logCat,'log_Iat'=dat$logIat,
                        'Mat'=Mat,'daysprop'=daysprop,
                        'tc'=tc.ssh,'robcode'=robcode)
parlist.robu.am <- c(opt.ml.am$par,list('log_Fat'=logFat.ini,
                                        'log_Nat'=logNat.ini))
obj.robu.am <- MakeADFun(data=datalist.rob.am,
                         parameters=parlist.robu.am,
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_st",silent=T)
system.time(opt.robu.am <- nlminb(start=obj.robu.am$par,obj=obj.robu.am$fn,
                                  gr=obj.robu.am$gr,
                                  control=list(eval.max=1000,iter.max=1000)))
# ^ about 2s on laptop
opt.robu.am
# ^ seems to have converged

system.time(rep.robu.am <- sdreport(obj.robu.am,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robu.am) # overwrite summ.rep, big object
logFat.robu.am <- matrix(summ.rep[which.logFat,1],AF,TF)
logNat.robu.am <- matrix(summ.rep[which.logNat,1],AN,TN)
# ^ predicted randeff based on uncorrected robust estimate, only for initial
#   values below


# Step 2/3 and 3/3: evaluate corrected score based on initial uncorrected robust
#                   estimate and correct it by single NR step

system.time(opt.robc.am <- NP_st_nrcorrect(theta.t=opt.robu.am$par,
                                            obj=obj.robu.am,H=H,maxit=1))
# ^ XXX seconds with H=1000 # TODO
str(opt.robc.am,1)
# ^ list with corrected robust gradient and (finite diff numerical) Hessian


# predict randeff based on corrected theta
parlist.robc <- c(opt.robc.am$theta.t,list('log_Fat'=logFat.robu.am,
                                           'log_Nat'=logNat.robu.am))
obj.robc.am <- MakeADFun(data=datalist.rob.am, # same data
                         parameters=parlist.robc, # updated par
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_st",silent=T)
invisible(obj.robc.am$fn())
# ^ no outer optim, just predict randeff with inner optim (Laplace approx)

system.time(rep.robc.am <- sdreport(obj.robc.am,bias.correct=F))
# ^ about 2 s on laptop

summ.rep <- summary(rep.robc.am) # overwrite summ.rep, big object
theta.robc.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1] # opt.robc.am$theta
theta.se.robc.am <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.robc.am <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.robc.am <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.robc.am <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.robc.am <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.robc.am <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.robc.am <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

w.rob.am <- NP_st_weights(obj.robc.am) # robustness weights
w.rob.am <- w.rob.am[5:8] # no need to rescale
str(w.rob.am)
# ^ list of 4 matrices, one for each loglik contribution: F_{a,t}, N_{a,t},
#   C_{a,t}, and I_{a,t}


# compare ML and robust estimates at the model
cbind(theta,theta.ml.am,theta.se.ml.am,theta.robc.am,theta.se.robc.am)
# ^ true theta, MLE and its se, robust and its se

par(mfrow=c(5,3))
for (j in 1:p){
  plot(1:2,c(theta.ml.am[j],theta.robc.am[j]),pch=19,
       xlim=c(0.7,2.3),ylim=theta[j]*c(0.5,1.5),xaxt='n',xlab='',ylab='',
       main=names.theta[j])
  arrows(x0=1:2,y0=c(theta.ml.am[j]-2*theta.se.ml.am[j],
                     theta.robc.am[j]-2*theta.se.robc.am[j]),
         x1=1:2,y1=c(theta.ml.am[j]+2*theta.se.ml.am[j],
                     theta.robc.am[j]+2*theta.se.robc.am[j]),
         angle=90,code=3,length=0.05)
  # ^ error bars as rough 95% CI, just to give a sense of variability
  axis(1,at=1:2,labels=c('ML','Robust'),lwd=0,lwd.ticks=1)
  abline(h=theta[j],col='red')
}
par(mfrow=c(1,1))
# ^ MLE and robust est generally nearly identical at the model


# compare ML and robust predicted randeff at the model

lb.Fat.ml.am <- exp(logFat.ml.am-2*logFat.se.ml.am)
ub.Fat.ml.am <- exp(logFat.ml.am+2*logFat.se.ml.am)
lb.Nat.ml.am <- exp(logNat.ml.am-2*logNat.se.ml.am)
ub.Nat.ml.am <- exp(logNat.ml.am+2*logNat.se.ml.am)
# ^ exponentiate bounds back on original scale

lb.Fat.robc.am <- exp(logFat.robc.am-2*logFat.se.robc.am)
ub.Fat.robc.am <- exp(logFat.robc.am+2*logFat.se.robc.am)
lb.Nat.robc.am <- exp(logNat.robc.am-2*logNat.se.robc.am)
ub.Nat.robc.am <- exp(logNat.robc.am+2*logNat.se.robc.am)
# ^ exponentiate bounds back on original scale

yearvec <- rep(NA,TF)
every5 <- (0:9)*5+4
yearvec[every5] <- year.F[every5] # label every 5 years

col.pollockages <- c('#da5e16','#7a13be','#f00000','#000cf0',
                     '#00850b','#c70074','#00c1d6','#6E7B8B')
col.pollockages.env <- paste0(col.pollockages,'30') # semi-transparent

names.age.F <- paste0('a = ',age.F)
names.age.N <- paste0('a = ',age.N)
names.age.I <- paste0('a = ',age.I)
names.age.C <- names.age.N



par(mfrow=c(2,2))
# st Fat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} ML predictions, at model',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.ml.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.ml.am[j,],ub.Fat.ml.am[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend('topright',names.age.F,lty=1,bty='n',
       col=col.pollockages[1:AF],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# st Fat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} robust predictions, at model',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.robc.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.robc.am[j,],ub.Fat.robc.am[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
# st Nat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} ML predictions, at model',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.ml.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.ml.am[j,],ub.Nat.ml.am[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend(x='topright',names.age.N,lty=1,bty='n',
       col=col.pollockages[1:AN],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# st Nat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} robust predictions, at model',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.robc.am[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.robc.am[j,],ub.Nat.robc.am[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
par(mfrow=c(1,1))


# investigate robustness weights at the model
lapply(w.rob.am,function(x){summary(as.numeric(x))})
# ^ nearly all weights are 1 as expected

unlist(lapply(w.rob.am,function(x){which(x<0.8,arr.ind=T)}))
# ^ only one really downweighted observation, given arbitrary threshold 0.8

dat$logCat[,11:15]
# ^ downweighted obs (6,13) may be atypically high given neighborhood


### contaminate simulated response
logCatcont <- dat$logCat
logCatcont[1,whichcontamC] <- pmax(logCatcont[1,whichcontamC]+addC,1)
# ^ contamination: underreport catch for recruits in year 2000, with floor at 1
#   on log scale to keep realistic values


### ML estimation, under contam
datalist.ml.uc <- list('log_Cat'=logCatcont,'log_Iat'=dat$logIat, # contaminated
                       'Mat'=Mat,'daysprop'=daysprop,
                       'tc'=rep(0,6),'robcode'=0)
# ^ data list fed to TMB's MakeADFun

obj.ml.uc <- MakeADFun(data=datalist.ml.uc, # contam data
                       parameters=parlist.ini, # same as at model
                       random=c('log_Nat','log_Fat'),
                       DLL="NP_st",silent=T)

# MLE: minimize (Laplace-approximated) marginal negloglik
system.time(opt.ml.uc <- nlminb(start=obj.ml.uc$par,obj=obj.ml.uc$fn,gr=obj.ml.uc$gr,
                                control=list(eval.max=1000,iter.max=1000)))
# ^ about 2s
opt.ml.uc
# ^ seems to have converged

# report theta estimates and randeff on original scale
system.time(rep.ml.uc <- sdreport(obj.ml.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.ml.uc)
theta.ml.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1]
theta.se.ml.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.ml.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.ml.uc <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.ml.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.ml.uc <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.ml.uc <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.ml.uc <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

cbind(theta,theta.ml.am,theta.ml.uc)
# ^ true theta, MLE at model and MLE under contam: sigmaF3-4 most biased


### Robust estimation, under contam

# Step 1/3: uncorrected robust estimator, minimize robustified marginal negloglik
datalist.rob.uc <- list('log_Cat'=logCatcont,'log_Iat'=dat$logIat, # contam
                        'Mat'=Mat,'daysprop'=daysprop,
                        'tc'=tc.ssh,'robcode'=robcode)
parlist.robu.uc <- c(opt.ml.uc$par,list('log_Fat'=logFat.ini,
                                        'log_Nat'=logNat.ini))
obj.robu.uc <- MakeADFun(data=datalist.rob.uc,
                         parameters=parlist.robu.uc,
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_st",silent=T)
system.time(opt.robu.uc <- nlminb(start=obj.robu.uc$par,obj=obj.robu.uc$fn,
                                  gr=obj.robu.uc$gr,
                                  control=list(eval.max=1000,iter.max=1000)))
# ^ about 3s
opt.robu.uc
# ^ seems to have converged

system.time(rep.robu.uc <- sdreport(obj.robu.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robu.uc) # overwrite summ.rep, big object
logFat.robu.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logNat.robu.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
# ^ predicted randeff based on uncorrected robust estimate, only for initial
#   values below


# Step 2/3 and 3/3: evaluate corrected score based on initial uncorrected robust
#                   estimate and correct it by single NR step

system.time(opt.robc.uc <- NP_st_nrcorrect(theta.t=opt.robu.uc$par,
                                            obj=obj.robu.uc,H=H,maxit=1))
# ^ XXX seconds with H=1000 # TODO
str(opt.robc.uc,1)
# ^ list with corrected robust gradient and (finite diff numerical) Hessian


# predict randeff based on corrected theta
parlist.robc <- c(opt.robc.uc$theta.t,list('log_Fat'=logFat.robu.uc,
                                           'log_Nat'=logNat.robu.uc))
obj.robc.uc <- MakeADFun(data=datalist.rob.uc, # same data
                         parameters=parlist.robc, # updated par
                         random=c('log_Nat','log_Fat'),
                         DLL="NP_st",silent=T)
invisible(obj.robc.uc$fn())
# ^ no outer optim, just predict randeff with inner optim (Laplace approx)

system.time(rep.robc.uc <- sdreport(obj.robc.uc,bias.correct=F))
# ^ about 2s

summ.rep <- summary(rep.robc.uc) # overwrite summ.rep, big object
theta.robc.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),1] # opt.robc.uc$theta
theta.se.robc.uc <- summ.rep[(p+AN*TN+AF*TF+1):(p+AN*TN+AF*TF+p),2]
logFat.robc.uc <- matrix(summ.rep[which.logFat,1],AF,TF)
logFat.se.robc.uc <- matrix(summ.rep[which.logFat,2],AF,TF)
logNat.robc.uc <- matrix(summ.rep[which.logNat,1],AN,TN)
logNat.se.robc.uc <- matrix(summ.rep[which.logNat,2],AN,TN)
Fat.robc.uc <- matrix(summ.rep[which.Fat,1],AF,TF) # pred on original scale
Nat.robc.uc <- matrix(summ.rep[which.Nat,1],AN,TN) # pred on original scale

w.rob.uc <- NP_st_weights(obj.robc.uc) # robustness weights
w.rob.uc <- w.rob.uc[5:8] # no need to rescale
str(w.rob.uc)
# ^ list of 4 matrices, one for each loglik contribution: F_{a,t}, N_{a,t},
#   C_{a,t}, and I_{a,t}


# compare ML and robust estimates under contam
cbind(theta,theta.ml.uc,theta.se.ml.uc,theta.robc.uc,theta.se.robc.uc)
# ^ true theta, MLE and its se, robust and its se

cbind(abs(theta.ml.am-theta.ml.uc)/theta,abs(theta.robc.am-theta.robc.uc)/theta)
# ^ est diff between at model and under contam, relative to true theta: MLE
#   changed a lot for a few param, robust remains globally stable

par(mfrow=c(5,3))
for (j in 1:p){
  plot(1:2,c(theta.ml.uc[j],theta.robc.uc[j]),pch=19,
       xlim=c(0.7,2.3),ylim=theta[j]*c(0.5,1.5),xaxt='n',xlab='',ylab='',
       main=names.theta[j])
  arrows(x0=1:2,y0=c(theta.ml.uc[j]-2*theta.se.ml.uc[j],
                     theta.robc.uc[j]-2*theta.se.robc.uc[j]),
         x1=1:2,y1=c(theta.ml.uc[j]+2*theta.se.ml.uc[j],
                     theta.robc.uc[j]+2*theta.se.robc.uc[j]),
         angle=90,code=3,length=0.05)
  # ^ error bars as rough 95% CI, just to give a sense of variability
  axis(1,at=1:2,labels=c('ML','Robust'),lwd=0,lwd.ticks=1)
  abline(h=theta[j],col='red')
}
par(mfrow=c(1,1))
# ^ robust est remains close to true value, MLE biased now


# compare ML and robust predicted randeff at the model
lb.Fat.ml.uc <- exp(logFat.ml.uc-2*logFat.se.ml.uc)
ub.Fat.ml.uc <- exp(logFat.ml.uc+2*logFat.se.ml.uc)
lb.Nat.ml.uc <- exp(logNat.ml.uc-2*logNat.se.ml.uc)
ub.Nat.ml.uc <- exp(logNat.ml.uc+2*logNat.se.ml.uc)
# ^ exponentiate bounds back on original scale

lb.Fat.robc.uc <- exp(logFat.robc.uc-2*logFat.se.robc.uc)
ub.Fat.robc.uc <- exp(logFat.robc.uc+2*logFat.se.robc.uc)
lb.Nat.robc.uc <- exp(logNat.robc.uc-2*logNat.se.robc.uc)
ub.Nat.robc.uc <- exp(logNat.robc.uc+2*logNat.se.robc.uc)
# ^ exponentiate bounds back on original scale


par(mfrow=c(2,2))
# st Fat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} ML predictions, under contam',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.ml.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.ml.uc[j,],ub.Fat.ml.uc[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend('topright',names.age.F,lty=1,bty='n',
       col=col.pollockages[1:AF],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# st Fat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(0,1.2),
     main='F_{a,t} robust predictions, under contam',
     ylab=expression('Fishing mortality'~italic(F)))
axis(side=1,at=1:TF,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TF,y=-0.1,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AF){
  lines(1:TF,Fat.robc.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TF,TF:1),
          y=c(lb.Fat.robc.uc[j,],ub.Fat.robc.uc[j,TF:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
# st Nat ML
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} ML predictions, under contam',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.ml.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.ml.uc[j,],ub.Nat.ml.uc[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
legend(x='topright',names.age.N,lty=1,bty='n',
       col=col.pollockages[1:AN],fill=paste0(col.pollockages[1:AF],10),
       border='transparent')
# st Nat Robust
plot(1:TF,type='n',xlab='',xaxt='n',ylim=c(600,2e6),
     main='N_{a,t} robust predictions, under contam',
     ylab=expression('Abundance'~italic(N)))
axis(side=1,at=1:TN,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1)
axis(side=1,at=every5,labels=NA,tick=T,lwd=0,line=0.1,lwd.ticks=1.5,tck=-0.03)
text(x=1:TN,y=-2e5,labels=yearvec,srt=45,adj=c(1,1),xpd=T)
for (j in 1:AN){
  lines(1:TN,Nat.robc.uc[j,],col=col.pollockages[j])
  polygon(x=c(1:TN,TN:1),
          y=c(lb.Nat.robc.uc[j,],ub.Nat.robc.uc[j,TN:1]),
          col=paste0(col.pollockages[j],10),border=NA,xpd=F)
}
par(mfrow=c(1,1))
# ^ ML predictions now biased around contminated year 2000 and trajectory of
#   recruits very jagged overall, robust predictions remain stable as expected


# investigate robustness weights under contamintion
lapply(w.rob.uc,function(x){summary(as.numeric(x))})
# ^ most weights are 1, but some very close to 0

lapply(w.rob.uc,function(x){which(x<0.8,arr.ind=T)})
# ^ one additional heavily downweighted commercial catch observation

logCatcont[,32:36]
# ^ downweighted obs (1,34) is exactly contaminated point in year 2000
dat$logCat[,32:36] # (1,34) was originally around 8























### END RobSSM_Main
