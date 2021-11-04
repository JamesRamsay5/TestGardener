Analyze <- function(theta, thetaQnt, dataList, ncycle=10, itdisp=FALSE, verbose=FALSE) {
  
  # Last modified 29 October 2021 by Jim Ramsay

  parList   <- vector("list",ncycle)  
  meanHsave <- rep(0,ncycle)
  
  # logdensbasis <- dataList$WfdPar$fd$basis
  logdensbasis <- create.bspline.basis(c(0,100), 15)
  
  WfdList <- dataList$WfdList
  
  parList <- vector("list", ncycle)
  
  U <- dataList$U
  
  for (icycle in 1:ncycle) {
    
    print(paste('Cycle ',icycle))
    
    #  ----------------------------------------------------------
    #  Step 1:  Bin the data, and smooth the binned data
    #  ----------------------------------------------------------
    
    if (verbose) print("Wbinsmth:")
    
    WfdResult <- Wbinsmth(theta, dataList, WfdList, thetaQnt)
    WfdList   <- WfdResult$WfdList
    binctr    <- WfdResult$aves
    bdry      <- WfdResult$bdry
    freq      <- WfdResult$freq
    
    #  compute current mean value of objective function H
    
    if (verbose) print("Hfun:")
    
    H <- Hfun(theta, WfdList, U)
    meanH <- mean(H)
    meanHsave[icycle] <- meanH
    print(paste('Mean surprisal = ', round(meanH,3)))
    
    #  ----------------------------------------------------------
    #  Step 2:  Compute optimal score index values
    #  ----------------------------------------------------------
    
    if (verbose) print("thetafun:")

    thetafunList <- thetafun(theta, WfdList, U, 20, 1e-3, itdisp)
    theta    <- thetafunList$theta_out
    Hval     <- thetafunList$Hval
    DHval    <- thetafunList$DHval
    D2Hval   <- thetafunList$D2Hval
    active   <- thetafunList$active
    
    #  ----------------------------------------------------------
    #  Step 3:  Estimate the score density for score index values
    #  ----------------------------------------------------------
    
    if (verbose) print("theta.distn:")
    
    thetadens <- theta[0 < theta & theta < 100]
    theta.distnList <- theta.distn(thetadens, logdensbasis)
    
    logdensfd <- theta.distnList$logdensfd
    denscdf   <- as.numeric(theta.distnList$denscdf)
    C         <- theta.distnList$C
    densfine  <- theta.distnList$densfine
    
    denscdfi <- unique(denscdf)
    indfinei <- seq(0,100,len=length(denscdfi))
    Qvec     <- pracma::interp1(denscdfi, indfinei, dataList$PcntMarkers/100)
    nbin     <- dataList$nbin
    thetaQnt <- pracma::interp1(denscdfi, indfinei, seq(0,2*nbin,1)/(2*nbin))
    thetaQnt[2*nbin+1] <- 100
    
    #  ----------------------------------------------------------
    #  Step 4.  Compute arc length and its measures
    #  ----------------------------------------------------------
      
    if (verbose) print("theta2arclen:")
    
    theta2arclenList <- theta2arclen(theta, WfdList, dataList$Wdim);
    theta_al      <- theta2arclenList$theta_al 
    arclength     <- theta2arclenList$arclength 
    arclengthfine <- theta2arclenList$arclengthfine
    Qvec_al       <- theta2arclenList$Qvec_al
    print(paste('arclength in bits = ',round(arclength,1)))
    
    #  ----------------------------------------------------------
    #  Step 5:  set up ParameterCell arrays
    #  ----------------------------------------------------------
    
    if (verbose) print("parList:")
    
    parListi <- list(
      theta      = theta,
      thetaQnt   = thetaQnt,
      WfdList    = WfdList,
      logdensfd  = logdensfd,
      C          = C,
      densfine   = densfine,
      denscdf    = denscdf,
      Qvec       = Qvec,
      binctr     = binctr,
      bdry       = bdry,
      freq       = freq,
      Hval       = Hval,
      DHval      = DHval,
      D2Hval     = D2Hval,
      active     = active,
      arclength  = arclength,
      alfine     = arclengthfine,
      Qvec_al    = Qvec_al,
      theta_al   = theta_al
    )
    
    parList[[icycle]] <- parListi
    
  }
  
  return(list(parList=parList, meanHsave=meanHsave))
  
} 