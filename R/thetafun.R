thetafun <- function(theta, WfdList, U, 
                     itermax = 20, crit = 1e-3, itdisp=FALSE) {
  
# Last modified 8 February 2021 by Jim Ramsay

  if (!inherits(WfdList,"list"))
    stop("Arguments WfdList is not list object.")
  
  N <- nrow(U)
  n <- ncol(U)
  
  if (length(theta) != N)
  {
    stop("Length of theta not equal to nrow(U).")
  }
  
  if (any(theta < 0) || any(theta > 100)) 
    stop('Values of theta out of bounds.')
  
  if (length(WfdList) != n)
    stop("Length of WfdList not equal to ncol(U).")
  
  # Initial function, first and second derivative values wrt theta
  
  Hval    <-  Hfun(theta, WfdList, U)
  DHList  <- DHfun(theta, WfdList, U)
  DHval   <- DHList$DH
  D2Hval  <- DHList$D2H
  
  # find better initial values for cases where D2H < 0
  # use minimum over a 26-value grid
  ngrid <- 26
  
  jneg <- which(D2Hval <= 0)
  if (length(jneg) > 0)
  {
    if (itdisp) {
      print(paste("Number of nonpositive D2H =",length(jneg)))
    }
    thetatry <- seq(0,100,length.out=ngrid)
    for (j in 1:length(jneg))
    {
      jind        <- jneg[j]
      Uj          <- matrix(U[jind,], 1)[rep(1,ngrid), ]
      Htry        <- Hfun(thetatry, WfdList, Uj)
      Hmin        <- min(Htry)
      kmin        <- which(Htry == Hmin)
      theta[jind] <- thetatry[kmin[1]]
    }
    Hvalnonpos   <-  Hfun(theta[jneg], WfdList, U[jneg,])
    result       <- DHfun(theta[jneg], WfdList, U[jneg,])
    DHvalnonpos  <- result$DH
    D2Hvalnonpos <- result$D2H
    Hval[jneg]   <- Hvalnonpos
    DHval[jneg]  <- DHvalnonpos
    D2Hval[jneg] <- D2Hvalnonpos
  }
  
  if (itdisp) {
    print(paste("mean adjusted values =",round(mean(theta[jneg]),2)))
  }
  
  if (itermax == 0)
  {
    theta_out <- theta
    return(theta_out)
  }
  
  # Initialize indices of active theta values
  
  active    <- ActiveTestFn(theta, DHval, D2Hval, crit)
  nactive   <- sum(active)
  
  # update scores in active set by a Newton-Raphson step
  
  thetaa <- theta[active]
  Ha     <- Hval[active]
  DHa    <- DHval[active]
  D2Ha   <- D2Hval[active]
  
  # loop through iterations to reduce number of active cases
  
  iter   <- 0
  while(nactive > 0 && iter < itermax)
  {
    iter <- iter + 1
    if (itdisp) {
      print(paste("iter",iter,",  nactive = ",nactive))
    }
    
    # compute Newton-Raphson step sizes
    step <- DHa/D2Ha
    
    # bound the steps between -10 and 10
    step[step < -10] <- -10
    step[step >  10] <-  10
    
    # ensure that initial step does not go to boundaries
    stepind0 <- which(thetaa - step <=   0)
    stepindn <- which(thetaa - step >= 100)
    if (length(stepind0) > 0)
    {
      step[stepind0] <-  (thetaa[stepind0] -   0)  
    }
    if (length(stepindn) > 0)
    {
      step[stepindn] <- -(100 - thetaa[stepindn])
    }
    
    #  half fac as required to achieve the reduction of all fn. values
    #  step through the iteration while any(Hval(active) - Ha < -1e-2)
    #  up to a maximum number of 10 step halvings
    thetaanew     <- theta[active]
    Uanew         <- as.matrix(U[active,])
    Hanew         <- Hval[active]
    fac           <- 1
    stepiter      <- 0
    
    #  initial set of theta's failing to meet criterion
    #  in this loop if initial step results in an increase in H,
    #  the step is halved until a reduction in H is achieved
    failindex  <- 1:nactive
    while(length(failindex) > 0 && stepiter < 10)
    {
      stepiter <- stepiter + 1
      #  compute current step size
      stepnew          <- fac * step[failindex]
      thetaanew[failindex]    <- thetaa[failindex] - stepnew
      #  Compute new function value
      if (nactive == 1)
      {
        Uanewfail <- Uanew
      } else
      {
        Uanewfail <- as.matrix(Uanew[failindex,])
      }
      # Compute new function value
      # print(failindex)
      # print(as.numeric(thetaanewfail))
      # print(which(is.na(thetaanewfail)))
      Hanew[failindex] <- Hfun(thetaanew[failindex], WfdList, Uanewfail) 
      failindex <- which((Hval[active] - Hanew) < -crit)
      #  halve fac in preparation for next step size iteration
      fac <- fac/2
    }
    
    # either all function values reduced or 10 halvings completed
    # save current function values
    
    theta[active] <- thetaanew
    
    # accept current value of H and also compute two derivatives
    Ha      <- Hanew
    DHLista <- DHfun(theta[active], WfdList, as.matrix(U[active,]))
    DHa     <- DHLista$DH
    D2Ha    <- DHLista$D2H
    jneg    <- which(D2Ha <= 0)
    
    Hval[active]   <- Ha
    DHval[active]  <- DHa
    D2Hval[active] <- D2Ha
    
    # find theta values that need further iterations
    # function ActiveTestFn has been modified
    active  <- ActiveTestFn(thetaanew, DHa, D2Ha, crit, active)
    nactive <- length(which(active == TRUE))
    if (nactive > 0)
    {
      thetaa  <- theta[active]
      DHa     <- DHval[active]
      D2Ha    <- D2Hval[active]
    }
    
    #  return for a new optimization step 
  }
  
  # opimization complete for all cases, save results and exit
  theta_out <- theta
  
  return(list(theta_out=theta_out, Hval=Hval, DHval=DHval, D2Hval=D2Hval, iter=iter))
  
}

#  ----------------------------------------------------------------------------

ActiveTestFn <- function(theta, DHval, D2Hval, crit = 1e-3, activeold = NULL) {
  # Last modified 29 January 2020 by Jim Ramsay
  
  if (is.null(activeold))  {
    N <- length(theta)
    activenew <- rep(FALSE, N)
    for (j in 1:N) {
      thetaj <- theta[j]
      test   <- (thetaj  >   0   && thetaj   < 100) |
                (thetaj ==   0   && DHval[j] <   0)   | 
                (thetaj == 100   && DHval[j] >   0)
      if (test) activenew[j] <- TRUE
    }
  } else  {
    nactive   <- length(theta)
    N         <- length(activeold)
    activenew <- rep(FALSE,N)
    activeind <- which(activeold)
    for (j in 1:nactive) {
      thetaj <- theta[j]
      DHaj   <- DHval[j]
      D2Haj  <- D2Hval[j]
      
      # 1:  abs(slope) > crit and D2H positive
      test1 <- thetaj > 10 && abs(DHaj) > 10 * crit && D2Haj > 0 && 
               thetaj >  0 && thetaj < 100
      
      # 2:  abs(slope) > crit and D2H positive
      test2 = thetaj <= 10 && abs(DHaj) > 50 * crit && D2Haj > 0 && 
              thetaj >   0 && thetaj < 100
      
      # 3:  theta at 100 and slope positive
      test3 <- thetaj == 0 && DHaj < 0
      
      # 4:  theta at 100 and slope negative
      test4 <- thetaj == 100 && DHaj > 0
      if (test1 || test2 || test3 || test4)
        activenew[activeind[j]] <- TRUE
    }
  }
  return(activenew)
}

