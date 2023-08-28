dataSimulation <- function(dataList, parList, theta.pop=seq(0,100,len=101), 
                           nsample=1000) {
  
  #  Simulate data, analyze each simulated sample, and print results
  
  #  Last modified on 12 November 2022
  
  #  --------------------------------------------------------------------------
  #                   Set up access to population info:
  #  --------------------------------------------------------------------------
  
  #  info from dataList
  
  ScoreList <- dataList$optList$optScr
  scrrng    <- dataList$scrrng
  Wdim      <- dataList$Wdim
  
  #  info from parList
  
  WfdList   <- parList$WfdList
  logdensfd <- parList$logdensfd
  Qvec      <- parList$Qvec
  binctr    <- parList$binctr
  n         <- length(WfdList)
  
  #  --------------------------------------------------------------------------
  #                       Generate simulated data
  #  --------------------------------------------------------------------------
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  #  define population theta values 
  
  ntheta  <- length(theta.pop)
  
  #  compute population values for theta, mu and arc length
  #  corresponding to nfine equally spaced values of theta.
  
  mu.pop        <- testscore(theta.pop, WfdList, dataList$optList)
  result        <- theta2arclen(theta.pop, Qvec, WfdList, binctr)
  al.pop        <- result$theta_al 
  arclength     <- result$arclength 
  arclength1001 <- result$arclength1001 
  arclengthfine <- seq(0,arclength,len=101)
  
  #  arrays to store simulated data
  
  Usave <- array(0,c(ntheta,n,nsample))
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) print(paste("Sample ",isample))
    Usave[,,isample] <- Usimulate(theta.pop, WfdList)
  }
  
  #  --------------------------------------------------------------------------
  #                    Analyze each sample
  #  --------------------------------------------------------------------------
  
  #  Note that total arc length is a function of the surprisal functions,
  #  and therefore is fixed for these simulations.
  
  #  set up matrices to save results
  
  sumscrsave <- matrix(0, ntheta, nsample)
  thetasave  <- matrix(0, ntheta, nsample)
  musave     <- matrix(0, ntheta, nsample)
  alsave     <- matrix(0, ntheta, nsample)
  
  Hvalsave   <- matrix(0, ntheta, nsample)
  DHvalsave  <- matrix(0, ntheta, nsample)
  D2Hvalsave <- matrix(0, ntheta, nsample)
  
  Infoopt <- FALSE  #  Use info matrix instead of Hessian
  
  #  analyze the samples
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) print(paste("Sample ", isample))
    # print(isample)
    Umati <- Usave[,,isample]
    #  sum scores for this sample
    # scrveci <- sumscorefn(Umati, ScoreList)
    scrveci <- matrix(0,ntheta,1)
    for (j in 1:ntheta) {
      for (i in 1:n) {
        scrveci[j] = scrveci[j] + ScoreList[[i]][Umati[j,i]]
      }
    }
    #  optimal scores for multi-option analysis
    result     <- thetafun(theta.pop, WfdList, Umati)             
    thetaveci  <- result$theta_out 
    Hvalveci   <- result$Hval 
    DHvalveci  <- result$DHval 
    D2Hvalveci <- result$D2Hval
    muveci     <- testscore(thetaveci, WfdList, dataList$optList)
    alveci     <- pracma::interp1(as.numeric(indfine), 
                                  as.numeric(arclengthfine), 
                                  as.numeric(thetaveci))
    sumscrsave[,isample] <-    scrveci
    thetasave[,isample]  <-  thetaveci     
    musave[,isample]     <-     muveci
    alsave[,isample]     <-     alveci
    Hvalsave[,isample]   <-   Hvalveci
    DHvalsave[,isample]  <-  DHvalveci
    D2Hvalsave[,isample] <- D2Hvalveci
  }
  
  #  --------------------------------------------------------------------------
  #             Define list  object simList to save the results
  #  --------------------------------------------------------------------------
  
  simList <- list(
    sumscr  = sumscrsave,
    theta   = thetasave,
    mu      = musave,
    al      = alsave,
    thepop  = theta.pop,
    mupop   = mu.pop,
    alpop   = al.pop,
    n       = n,
    ntheta  = ntheta,
    indfine = indfine, 
    Qvec    = Qvec
  )
  
  return(simList)
  
}
