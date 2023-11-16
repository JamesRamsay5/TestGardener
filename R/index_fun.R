index_fun <- function(index, SfdList, chcemat, 
                     itermax = 20, crit = 1e-3, itdisp=FALSE) {
  
# Last modified 2 November 23 by Jim Ramsay

  if (!inherits(SfdList,"list"))
    stop("Arguments SfdList is not list object.")
  
  N <- nrow(chcemat)
  n <- ncol(chcemat)
  
  if (length(index) != N)
  {
    stop("Length of index not equal to nrow(chcemat).")
  }
  
  if (any(index < 0) || any(index > 100)) 
    stop('Values of index out of bounds.')
  
  if (length(SfdList) != n)
    stop("Length of SfdList not equal to ncol(chcemat).")
  
  # Initial function, first and second derivative values wrt index
  
  Fval    <-  Ffun(index, SfdList, chcemat)
  DFList  <- DFfun(index, SfdList, chcemat)
  DFval   <- DFList$DF
  D2Fval  <- DFList$D2F
  
  # find better initial values for cases where D2F < 0
  # use minimum over a 26-value grid
  ngrid <- 26
  
  jneg <- which(D2Fval <= 0)
  if (length(jneg) > 0)
  {
    if (itdisp) {
      print(paste("Number of nonpositive D2F =",length(jneg)))
    }
    indextry <- seq(0,100,length.out=ngrid)
    for (j in 1:length(jneg))
    {
      jind        <- jneg[j]
      chcematj          <- matrix(chcemat[jind,], 1)[rep(1,ngrid), ]
      Ftry        <- Ffun(indextry, SfdList, chcematj)
      Fmin        <- min(Ftry)
      kmin        <- which(Ftry == Fmin)
      index[jind] <- indextry[kmin[1]]
    }
    Fvalnonpos   <-  Ffun(index[jneg], SfdList, chcemat[jneg,])
    result       <- DFfun(index[jneg], SfdList, chcemat[jneg,])
    DFvalnonpos  <- result$DF
    D2Fvalnonpos <- result$D2F
    Fval[jneg]   <- Fvalnonpos
    DFval[jneg]  <- DFvalnonpos
    D2Fval[jneg] <- D2Fvalnonpos
  }
  
  if (itdisp) {
    print(paste("mean adjusted values =",round(mean(index[jneg]),2)))
  }
  
  if (itermax == 0)
  {
    index_out <- index
    return(index_out)
  }
  
  # Initialize indices of active index values
  
  active    <- ActiveTestFn(index, DFval, D2Fval, crit)
  nactive   <- sum(active)
  
  # update scores in active set by a Newton-Raphson step
  
  indexa <- index[active]
  Fa     <- Fval[active]
  DFa    <- DFval[active]
  D2Fa   <- D2Fval[active]
  
  # loop through iterations to reduce number of active cases
  
  iter   <- 0
  while(nactive > 0 && iter < itermax)
  {
    iter <- iter + 1
    if (itdisp) {
      print(paste("iter",iter,",  nactive = ",nactive))
    }
    
    # compute Newton-Raphson step sizes
    step <- DFa/D2Fa
    
    # bound the steps between -10 and 10
    step[step < -10] <- -10
    step[step >  10] <-  10
    
    # ensure that initial step does not go to boundaries
    stepind0 <- which(indexa - step <=   0)
    stepindn <- which(indexa - step >= 100)
    if (length(stepind0) > 0)
    {
      step[stepind0] <-  (indexa[stepind0] -   0)  
    }
    if (length(stepindn) > 0)
    {
      step[stepindn] <- -(100 - indexa[stepindn])
    }
    
    #  half fac as required to achieve the reduction of all fn. values
    #  step through the iteration while any(Fval(active) - Fa < -1e-2)
    #  up to a maximum number of 10 step halvings
    indexanew     <- index[active]
    chcematanew         <- as.matrix(chcemat[active,])
    Fanew         <- Fval[active]
    fac           <- 1
    stepiter      <- 0
    
    #  initial set of index's failing to meet criterion
    #  in this loop if initial step results in an increase in F,
    #  the step is halved until a reduction in F is achieved
    failindex  <- 1:nactive
    while(length(failindex) > 0 && stepiter < 10)
    {
      stepiter <- stepiter + 1
      #  compute current step size
      stepnew          <- fac * step[failindex]
      indexanew[failindex]    <- indexa[failindex] - stepnew
      #  Compute new function value
      if (nactive == 1)
      {
        chcematanewfail <- chcematanew
      } else
      {
        chcematanewfail <- as.matrix(chcematanew[failindex,])
      }
      # Compute new function value
      # print(failindex)
      # print(as.numeric(indexanewfail))
      # print(which(is.na(indexanewfail)))
      Fanew[failindex] <- Ffun(indexanew[failindex], SfdList, chcematanewfail) 
      failindex <- which((Fval[active] - Fanew) < -crit)
      #  halve fac in preparation for next step size iteration
      fac <- fac/2
    }
    
    # either all function values reduced or 10 halvings completed
    # save current function values
    
    index[active] <- indexanew
    
    # accept current value of F and also compute two derivatives
    Fa      <- Fanew
    DFLista <- DFfun(index[active], SfdList, as.matrix(chcemat[active,]))
    DFa     <- DFLista$DF
    D2Fa    <- DFLista$D2F
    jneg    <- which(D2Fa <= 0)
    
    Fval[active]   <- Fa
    DFval[active]  <- DFa
    D2Fval[active] <- D2Fa
    
    # find index values that need further iterations
    # function ActiveTestFn has been modified
    active  <- ActiveTestFn(indexanew, DFa, D2Fa, crit, active)
    nactive <- length(which(active == TRUE))
    if (nactive > 0)
    {
      indexa  <- index[active]
      DFa     <- DFval[active]
      D2Fa    <- D2Fval[active]
    }
    
    #  return for a new optimization step 
  }
  
  # opimization complete for all cases, save results and exit
  index_out <- index
  
  return(list(index_out=index_out, Fval=Fval, DFval=DFval, D2Fval=D2Fval, iter=iter))
  
}

#  ----------------------------------------------------------------------------

ActiveTestFn <- function(index, DFval, D2Fval, crit = 1e-3, activeold = NULL) {
  # Last modified 29 January 2020 by Jim Ramsay
  
  if (is.null(activeold))  {
    N <- length(index)
    activenew <- rep(FALSE, N)
    for (j in 1:N) {
      indexj <- index[j]
      test   <- (indexj  >   0   && indexj   < 100) |
                (indexj ==   0   && DFval[j] <   0)   | 
                (indexj == 100   && DFval[j] >   0)
      if (test) activenew[j] <- TRUE
    }
  } else  {
    nactive   <- length(index)
    N         <- length(activeold)
    activenew <- rep(FALSE,N)
    activeind <- which(activeold)
    for (j in 1:nactive) {
      indexj <- index[j]
      DFaj   <- DFval[j]
      D2Faj  <- D2Fval[j]
      
      # 1:  abs(slope) > crit and D2F positive
      test1 <- indexj > 10 && abs(DFaj) > 10 * crit && D2Faj > 0 && 
               indexj >  0 && indexj < 100
      
      # 2:  abs(slope) > crit and D2F positive
      test2 = indexj <= 10 && abs(DFaj) > 50 * crit && D2Faj > 0 && 
              indexj >   0 && indexj < 100
      
      # 3:  index at 100 and slope positive
      test3 <- indexj == 0 && DFaj < 0
      
      # 4:  index at 100 and slope negative
      test4 <- indexj == 100 && DFaj > 0
      if (test1 || test2 || test3 || test4)
        activenew[activeind[j]] <- TRUE
    }
  }
  return(activenew)
}

