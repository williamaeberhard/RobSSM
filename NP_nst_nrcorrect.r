NP_nst_nrcorrect <- function(theta.t=NULL,obj,H=100,seedval=NULL,
                             update.hessian=F,correct.hessian=F,
                             onlyimprove=T,
                             trace=T,maxit=1,tol.gr=1e-3){
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Setup ####
  #/////////////////////////////////////////////////////////////////////////////
  # adapted from TMB::checkConsistency
  
  args <- as.list(obj$env)[intersect(names(formals(MakeADFun)),ls(obj$env))]
  r0 <- r <- obj$env$random
  if (is.null(theta.t)){
    if (is.null(obj$env$last.par.best)){
      stop('"theta.t" not specified.')
    }
    parfull <- obj$env$last.par.best # if not supplied, take best param
    if (any(r)){
      theta.t <- parfull[-r] # retain only fixed param
    } else {
      theta.t <- parfull # if no randeff
    }
  } else {
    if (is.null(obj$env$last.par.best)){
      parfull <- obj$env$par # initial values
    } else {
      parfull <- obj$env$last.par.best # take best param
    }
    if (any(r)){
      parfull[-r] <- theta.t # replace fixed param by supplied value
    } else {
      parfull <- theta.t # if no randeff
    }
  }
  if (any(obj$env$profile)){
    r0 <- r[!as.logical(obj$env$profile)]
    names.profile <- unique(names(parfull[r[as.logical(obj$env$profile)]]))
  } else {
    names.profile <- NULL
  }
  names.random <- unique(names(parfull[r0]))
  args$parameters <- obj$env$parList(theta.t, par=parfull)
  # ^ replaces ini by theta.t and best randeff
  map.profile <- lapply(args$parameters[names.profile],function(x){factor(x*NA)})
  args$map <- c(args$map, map.profile)
  args$random <- names.random
  args$regexp <- FALSE
  newobj <- do.call("MakeADFun", args)
  # ^ create newobj similar to supplied obj, with theta.t if supplied
  obj.h <- do.call("MakeADFun", args) # similar, but only for simulation
  
  which.Fat <- which(names(newobj$env$par)=='log_Fat')
  which.Nat <- which(names(newobj$env$par)=='log_Nat')
  AF <- nrow(newobj$env$parameters$log_Fat) # 3, ..., 9+
  TF <- ncol(newobj$env$parameters$log_Fat) # t = 1967, ..., 2016
  AN <- nrow(newobj$env$parameters$log_Nat) # 3, ..., 10+
  TN <- ncol(newobj$env$parameters$log_Nat) # t = 1967, ..., 2016
  theta.t.upd <- theta.t # do no overwrite initial theta.t, to report back
  theta.upd <- NP_nst_untransfo(theta.t.upd)
  
  # extract predicted ini logNat and logFat on original data, fixed for simul
  logF1.fixed <- newobj$env$last.par.best[which.Fat][1:AF]
  logN1.fixed <- newobj$env$last.par.best[which.Nat][1:AN]
  
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Eval initial corrgr (and corrhess) ####
  #/////////////////////////////////////////////////////////////////////////////
  
  # eval uncorrected gr and he on original data, inital param
  gr.rob <- newobj$gr(theta.t.upd)/TF # sample size = time series length
  he.rob <- optimHess(par=theta.t.upd,fn=newobj$fn,gr=newobj$gr)/TF
  
  if (correct.hessian){
    # eval gr and he on simulated data => Egrad and Ehess
    gradmat <- matrix(NA_real_,H,length(theta.upd))
    Ehess <- matrix(0,length(theta.upd),length(theta.upd))
    H.ok <- 0
    if (!is.null(seedval)){set.seed(seedval)}
    for (h in 1:H){
      samlist.h <- NP_nst_gen(theta=theta.upd,
                              logF1=logF1.fixed,logN1=logN1.fixed)
      obj.h$env$data$log_Cat <- samlist.h$logCat
      obj.h$env$data$log_Iat <- samlist.h$logIat
      # ^ replace by simulated data
      
      obj.h$env$last.par.best[which.Fat] <- as.numeric(samlist.h$logFat)
      obj.h$env$last.par.best[which.Nat] <- as.numeric(samlist.h$logNat)
      # ^ replace randeff ini by true randeff, helps inner Newton
      
      obj.h$env$L.created.by.newton <- NULL # re do Cholesky
      obj.h$env$retape() # runs Newton for new pred randeff
      
      tryit1 <- try(obj.h$gr(theta.t.upd),T) # eval grad, simulated data
      if (any(is(tryit1,'try-error') | !is.finite(tryit1))){
        # if (trace){warning('h = ',h,' failed on gr')} # rarely fails anyway
      } else {
        tryit2 <- try(optimHess(par=theta.t.upd,fn=obj.h$fn,gr=obj.h$gr),T)
        if (any(is(tryit2,'try-error') | !is.finite(tryit2))){
          # if (trace){warning('h = ',h,' failed on he')}
        } else {
          gradmat[h,] <- tryit1
          Ehess <- Ehess+tryit2
          H.ok <- H.ok+1
        }
      }
    }
    Egrad <- colSums(na.omit(gradmat))/H.ok/TF # empirical average of gr
    corrgr <- gr.rob-Egrad # corrected robust gradient
    Ehess <- Ehess/H.ok/TF # empirical average of he
    corrhess <- he.rob-Ehess # corrected robust Hessian (1st order approx)
  } else {
    # eval only gr on simulated data => Egrad
    gradmat <- matrix(NA_real_,H,length(theta.upd))
    H.ok <- 0
    if (!is.null(seedval)){set.seed(seedval)}
    for (h in 1:H){
      samlist.h <- NP_nst_gen(theta=theta.upd,
                              logF1=logF1.fixed,logN1=logN1.fixed)
      obj.h$env$data$log_Cat <- samlist.h$logCat
      obj.h$env$data$log_Iat <- samlist.h$logIat
      # ^ replace by simulated data
      
      obj.h$env$last.par.best[which.Fat] <- as.numeric(samlist.h$logFat)
      obj.h$env$last.par.best[which.Nat] <- as.numeric(samlist.h$logNat)
      # ^ replace randeff ini by true randeff, helps inner Newton
      
      obj.h$env$L.created.by.newton <- NULL # re do Cholesky
      obj.h$env$retape() # runs Newton for new pred randeff
      
      tryit1 <- try(obj.h$gr(theta.t.upd),T) # eval grad, simulated data
      if (any(is(tryit1,'try-error') | !is.finite(tryit1))){
        # if (trace){warning('h = ',h,' failed on gr')} # rarely fails anyway
      } else {
        gradmat[h,] <- tryit1
        H.ok <- H.ok+1
      }
    }
    Egrad <- colSums(na.omit(gradmat))/H.ok/TF # empirical average of gr
    corrgr <- gr.rob-Egrad # corrected robust gradient
  }
  
  if (trace){
    message('it = 0: max(abs(corrgr)) = ',round(max(abs(corrgr)),6))
    cat('H.ok =',H.ok,'out of H =',H,'\n')
    cat('theta =',round(theta.upd,3),'\n')
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Iterate (pseudo) NR steps ####
  #/////////////////////////////////////////////////////////////////////////////
  
  if (update.hessian){
    if (correct.hessian){
      # update and correct Hessian at each it
      it <- 1
      while (max(abs(corrgr))>=tol.gr & it<=maxit){
        
        # keep previous step in case corrgr worsens, neglect Hessian for now
        corrgr0 <- corrgr
        theta.t.upd0 <- theta.t.upd
        
        # correct theta.t by NR step
        invcorrhess <- solve(corrhess)
        theta.t.upd <- theta.t.upd-as.numeric(invcorrhess%*%corrgr) # NR
        theta.upd <- NP_nst_untransfo(theta.t.upd)
        
        # update newobj and obj.h with new param
        newobj$env$parameters[-r] <- as.list(theta.t.upd)
        newobj$env$par[-r] <- theta.t.upd
        newobj$env$par.last.best[-r] <- theta.t.upd
        obj.h$env$parameters[-r] <- as.list(theta.t.upd)
        obj.h$env$par[-r] <- theta.t.upd
        obj.h$env$par.last.best[-r] <- theta.t.upd
        # ^ are all three necessary?
        
        # update uncorr gr and he on original data
        newobj$env$L.created.by.newton <- NULL # redo Cholesky
        newobj$env$retape() # not sure if necessary here
        gr.rob <- newobj$gr(theta.t.upd)/TF # sample size = time series length
        
        tryit.he <- try(optimHess(par=theta.t.upd,fn=newobj$fn,
                                  gr=newobj$gr)/TF,T)
        if (all(!is(tryit.he,'try-error'))){
          he.rob <- tryit.he
        } # else keep previous he.rob
        
        # eval gr and he on simulated data => Egrad and Ehess
        gradmat <- matrix(NA_real_,H,length(theta.upd))
        Ehess <- matrix(0,length(theta.upd),length(theta.upd))
        H.ok <- 0
        if (!is.null(seedval)){set.seed(seedval)}
        for (h in 1:H){
          samlist.h <- NP_nst_gen(theta=theta.upd,
                                  logF1=logF1.fixed,logN1=logN1.fixed)
          obj.h$env$data$log_Cat <- samlist.h$logCat
          obj.h$env$data$log_Iat <- samlist.h$logIat
          # ^ replace by simulated data
          
          obj.h$env$last.par.best[which.Nat] <- as.numeric(samlist.h$logNat)
          obj.h$env$last.par.best[which.Fat] <- as.numeric(samlist.h$logFat)
          # ^ replace randeff ini by true randeff, helps inner Newton
          
          obj.h$env$L.created.by.newton <- NULL # re do Cholesky
          obj.h$env$retape() # runs Newton for new pred randeff
          
          tryit1 <- try(obj.h$gr(theta.t.upd),T) #  eval grad only simulated data
          if (any(is(tryit1,'try-error') | !is.finite(tryit1))){
            # if (trace){warning('h = ',h,' failed on gr')} # rarely fails anyway
          } else {
            tryit2 <- try(optimHess(par=theta.t.upd,fn=obj.h$fn,gr=obj.h$gr),T)
            if (any(is(tryit2,'try-error') | !is.finite(tryit2))){
              # if (trace){warning('h = ',h,' failed on he')}
            } else {
              gradmat[h,] <- tryit1
              Ehess <- Ehess+tryit2
              H.ok <- H.ok+1
            }
          }
        }
        Egrad <- colSums(na.omit(gradmat))/H.ok/TF # empirical average of gr
        corrgr <- gr.rob-Egrad # corrected robust gradient
        Ehess <- Ehess/H.ok/TF # empirical average of he
        corrhess <- he.rob-Ehess # corrected robust Hessian (1st order approx)
        
        if (onlyimprove){
          which.worse <- which(abs(corrgr0)<abs(corrgr)) # NR worsened corrgr
          if (length(which.worse)!=0){ # otherwise keep theta.upd as is
            theta.t.upd[which.worse] <- theta.t.upd0[which.worse] # keep previous
            theta.upd <- NP_nst_untransfo(theta.t.upd)
            corrgr[which.worse] <- corrgr0[which.worse]
          }
        }
        
        if (trace){
          message('it = ',it,': max(abs(corrgr)) = ',round(max(abs(corrgr)),6))
          cat('H.ok =',H.ok,'out of H =',H,'\n')
          cat('theta =',round(theta.upd,3),'\n')
        }
        
        it <- it+1
      }
    } else {
      # update Hessian at each it, but do not correct it
      it <- 1
      while (max(abs(corrgr))>=tol.gr & it<=maxit){
        
        # keep previous step in case corrgr worsens
        corrgr0 <- corrgr
        theta.t.upd0 <- theta.t.upd
        
        # correct theta.t by pseudo-NR step
        invcorrhess <- solve(he.rob)
        theta.t.upd <- theta.t.upd-as.numeric(invcorrhess%*%corrgr) # NR
        theta.upd <- NP_nst_untransfo(theta.t.upd)
        
        # update newobj and obj.h with new param
        newobj$env$parameters[-r] <- as.list(theta.t.upd)
        newobj$env$par[-r] <- theta.t.upd
        newobj$env$par.last.best[-r] <- theta.t.upd
        obj.h$env$parameters[-r] <- as.list(theta.t.upd)
        obj.h$env$par[-r] <- theta.t.upd
        obj.h$env$par.last.best[-r] <- theta.t.upd
        # ^ are all three necessary?
        
        # update uncorr gr and he on original data
        newobj$env$L.created.by.newton <- NULL # redo Cholesky
        newobj$env$retape() # not sure if necessary here
        gr.rob <- newobj$gr()/TF # sample size = time series length
        tryit.he <- try(optimHess(par=theta.t.upd,fn=newobj$fn,
                                  gr=newobj$gr)/TF,T)
        if (all(!is(tryit.he,'try-error'))){
          he.rob <- tryit.he
        } # else keep previous he.rob
        
        # eval gr on simulated data => Egrad
        gradmat <- matrix(NA_real_,H,length(theta.upd))
        H.ok <- 0
        if (!is.null(seedval)){set.seed(seedval)}
        for (h in 1:H){
          samlist.h <- NP_nst_gen(theta=theta.upd,
                                  logF1=logF1.fixed,logN1=logN1.fixed)
          obj.h$env$data$log_Cat <- samlist.h$logCat
          obj.h$env$data$log_Iat <- samlist.h$logIat
          # ^ replace by simulated data
          
          obj.h$env$last.par.best[which.Fat] <- as.numeric(samlist.h$logFat)
          obj.h$env$last.par.best[which.Nat] <- as.numeric(samlist.h$logNat)
          # ^ replace randeff ini by true randeff, helps inner Newton
          
          obj.h$env$L.created.by.newton <- NULL # re do Cholesky
          obj.h$env$retape() # runs Newton for new pred randeff
          
          tryit1 <- try(obj.h$gr(theta.t.upd),T) #  eval grad only simulated data
          if (any(is(tryit1,'try-error') | !is.finite(tryit1))){
            # if (trace){warning('h = ',h,' failed on gr')} # rarely fails anyway
          } else {
            gradmat[h,] <- tryit1
            H.ok <- H.ok+1
          }
        }
        Egrad <- colSums(na.omit(gradmat))/H.ok/TF # empirical average of gr
        corrgr <- gr.rob-Egrad # corrected robust gradient
        
        if (onlyimprove){
          which.worse <- which(abs(corrgr0)<abs(corrgr)) # NR worsened corrgr
          if (length(which.worse)!=0){ # otherwise keep theta.upd as is
            theta.t.upd[which.worse] <- theta.t.upd0[which.worse] # keep previous
            theta.upd <- NP_nst_untransfo(theta.t.upd)
            corrgr[which.worse] <- corrgr0[which.worse]
          }
        }
        
        if (trace){
          message('it = ',it,': max(abs(corrgr)) = ',round(max(abs(corrgr)),6))
          cat('H.ok =',H.ok,'out of H =',H,'\n')
          cat('theta =',round(theta.upd,3),'\n')
        }
        it <- it+1
      }
    }
  } else { # do not update Hessian
    if (correct.hessian){
      # keep initial corrected hessian
      invcorrhess <- solve(corrhess)
    } else {
      # keep initial uncorrected hessian
      invcorrhess <- solve(he.rob)
    }
    
    it <- 1
    while (max(abs(corrgr))>=tol.gr & it<=maxit){
      
      # keep previous step in case corrgr worsens
      corrgr0 <- corrgr
      theta.t.upd0 <- theta.t.upd
      
      if (onlyimprove){
        # don't update if grad close to 0, unstable Hessian might mess up
        stepsize <- as.numeric(!(abs(corrgr)<tol.gr &
                                   apply(abs(gradmat),2,mean,na.rm=T,trim=0)<1))
        # correct theta.t by pseudo-NR step, only elements with grad far from 0
        theta.t.upd <- theta.t.upd-stepsize*as.numeric(invcorrhess%*%corrgr) # NR
      } else {
        # correct theta.t by pseudo-NR step, all elements
        theta.t.upd <- theta.t.upd-as.numeric(invcorrhess%*%corrgr) # NR
      }
      theta.upd <- NP_nst_untransfo(theta.t.upd)
      
      # update newobj and obj.h with new param
      newobj$env$parameters[-r] <- as.list(theta.t.upd)
      newobj$env$par[-r] <- theta.t.upd
      newobj$env$par.last.best[-r] <- theta.t.upd
      obj.h$env$parameters[-r] <- as.list(theta.t.upd)
      obj.h$env$par[-r] <- theta.t.upd
      obj.h$env$par.last.best[-r] <- theta.t.upd
      # ^ are all three necessary?
      
      # update uncorr gr on original data
      newobj$env$L.created.by.newton <- NULL # not sure if necessary here
      newobj$env$retape() # not sure if necessary here
      gr.rob <- newobj$gr(theta.t.upd)/TF # sample size = time series length
      
      # eval gr on simulated data => Egrad
      gradmat <- matrix(NA_real_,H,length(theta.upd))
      H.ok <- 0
      if (!is.null(seedval)){set.seed(seedval)}
      for (h in 1:H){
        samlist.h <- NP_nst_gen(theta=theta.upd,
                                logF1=logF1.fixed,logN1=logN1.fixed)
        obj.h$env$data$log_Cat <- samlist.h$logCat
        obj.h$env$data$log_Iat <- samlist.h$logIat
        # ^ replace by simulated data
        
        obj.h$env$last.par.best[which.Fat] <- as.numeric(samlist.h$logFat)
        obj.h$env$last.par.best[which.Nat] <- as.numeric(samlist.h$logNat)
        # ^ replace randeff ini by true randeff, helps inner Newton
        
        obj.h$env$L.created.by.newton <- NULL # re do Cholesky
        obj.h$env$retape() # runs Newton for new pred randeff
        
        tryit1 <- try(obj.h$gr(theta.t.upd),T) # eval grad, simulated data
        if (any(is(tryit1,'try-error') | !is.finite(tryit1))){
          # if (trace){warning('h = ',h,' failed on gr')} # rarely fails anyway
        } else {
          gradmat[h,] <- tryit1
          H.ok <- H.ok+1
        }
      }
      Egrad <- colSums(na.omit(gradmat))/H.ok/TF # empirical average of gr
      corrgr <- gr.rob-Egrad # corrected robust gradient
      
      if (onlyimprove){
        which.worse <- which(abs(corrgr0)<abs(corrgr)) # NR worsened corrgr
        if (length(which.worse)!=0){ # otherwise keep theta.upd as is
          theta.t.upd[which.worse] <- theta.t.upd0[which.worse] # keep previous
          theta.upd <- NP_nst_untransfo(theta.t.upd)
          corrgr[which.worse] <- corrgr0[which.worse]
        }
      }
      
      if (trace){
        message('it = ',it,': max(abs(corrgr)) = ',round(max(abs(corrgr)),6))
        cat('H.ok =',H.ok,'out of H =',H,'\n')
        cat('theta =',round(theta.upd,3),'\n')
      }
      it <- it+1
    }
  }
  
  #/////////////////////////////////////////////////////////////////////////////
  #### Output ####
  #/////////////////////////////////////////////////////////////////////////////
  
  if (!correct.hessian){
    corrhess <- NA
  }
  res <- list('theta'=theta.upd,'theta.t'=theta.t.upd,'theta.t.ini'=theta.t,
              'gradient.corrected'=corrgr,
              'hessian.corrected'=corrhess,'hessian.uncorrected'=he.rob)
  return(res)
}
# END NP_nst_nrcorrect
