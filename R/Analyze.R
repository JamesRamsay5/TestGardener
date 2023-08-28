Analyze <- function(theta, thetaQnt, dataList, ncycle=10, itdisp=FALSE, 
                    verbose=FALSE) {
  
  # Last modified 21 February 2023 by Jim Ramsay

  #  set up list vector to contain all results for each cycle
  
  parList <- vector("list",ncycle)  
  
  #  define the spline basis for representing the log density function
  
  logdensbasis <- fda::create.bspline.basis(c(0,100), 15)
  
  #  initialize surprisal curves and compute size of overspace
  
  WfdList <- dataList$WfdList
  
  n = length(WfdList)
  
  Wdim = 0
  for (i in 1:n) {
    WStri = WfdList[[i]]
    Wdim  = Wdim + WStri$M
  }
  
  #  set up data matrix
  
  U <- dataList$U
  
  #  main cycle loop
  
  for (icycle in 1:ncycle) {
    
    #  if (verbose) print(paste('----------  Cycle ',icycle,'-----------'))
    print(paste('----------  Cycle ',icycle,'-----------'))
    
    #  ----------------------------------------------------------
    #  Step 1:  Bin the data, and smooth the binned data
    #  ----------------------------------------------------------
    # print("step 1")
    if (verbose) print("Optimize surprisal curves:")
    
    WfdResult <- Wbinsmth(theta, dataList, WfdList, thetaQnt)
    WfdList   <- WfdResult$WfdList
    
    #  ----------------------------------------------------------
    #  Step 2:  compute mean value of objective function H
    #  ----------------------------------------------------------
    
    if (verbose) print("Compute mean examinee fits")
    
    H <- Hfun(theta, WfdList, U)
    meanH <- mean(H)
    
    if (verbose) print(paste('Mean data fit = ', round(meanH,3)))
    
    #  ----------------------------------------------------------
    #  Step 3:  Compute optimal score index values
    #  ----------------------------------------------------------
    # print("step 3")
    
    if (verbose) print("Optimize examinee data fits")

    thetafunList <- thetafun(theta, WfdList, U, 20, 1e-3, itdisp=itdisp)
    theta    <- thetafunList$theta_out
    Hval     <- thetafunList$Hval
    DHval    <- thetafunList$DHval
    D2Hval   <- thetafunList$D2Hval
    active   <- thetafunList$active
 
    #  ----------------------------------------------------------
    #  Step 4:  Estimate the score density for score index values
    #  ----------------------------------------------------------
    # print("step 4")
    
    if (verbose) print("Compute score index density")
    
    thetadens <- theta[0 < theta & theta < 100]
    theta.distnList <- theta.distn(thetadens, logdensbasis)
    pdf_fd    <- theta.distnList$pdf_fd
    logdensfd <- theta.distnList$logdensfd
    cdffine   <- theta.distnList$cdffine
    C         <- theta.distnList$C
    denscdf   <- theta.distnList$denscdf
    indcdf    <- theta.distnList$indcdf
    markers   <- dataList$PcntMarkers/100
    Qvec      <- pracma::interp1(as.numeric(denscdf), as.numeric(indcdf), markers)
    nbin      <- dataList$nbin
    bdry      <- seq(0,2*nbin,1)/(2*nbin)
    thetaQnt  <- pracma::interp1(as.numeric(denscdf), as.numeric(indcdf), bdry)
    
    #  ----------------------------------------------------------
    #  Step 5.  Compute arc length and its measures
    #  ----------------------------------------------------------
    
    DWfine = matrix(0,101,Wdim)
    m2 = 0
    for (i in 1:n) {
      WListi = WfdList[[i]]
      Mi     = WListi$M
      m1 = m2 + 1
      m2 = m2 + Mi
      DWfine[,m1:m2] = WListi$DWmatfine
    }
    indfine <- seq(0,100,len=101)
    arclength = max(pracma::cumtrapz(indfine, sqrt(apply(DWfine^2,1,sum))))
    
    if (verbose)  print(paste('arclength in bits = ',round(arclength,1)))
    
    #  ----------------------------------------------------------
    #  Step 6:  Check for mis-identications of minimum theta 
    #  ----------------------------------------------------------
    
    if (dataList$WfdPar$fd$basis$nbasis > 2) {
      Result <- thetasearch(WfdList, dataList$U, theta, Hval, DHval, D2Hval)
      theta  <- Result$theta
      Hval   <- Result$Hval
      DHval  <- Result$DHval
      D2Hval <- Result$D2Hval
    }
    print("thetasearch done")
    
    #  ----------------------------------------------------------
    #  Step 7:  set up ParameterCell arrays
    #  ----------------------------------------------------------

    parListi <- list(
      theta      = theta,
      thetaQnt   = thetaQnt,
      WfdList    = WfdResult$WfdList,
      meanH      = meanH,
      binctr     = WfdResult$aves,
      bdry       = WfdResult$bdry,
      freq       = WfdResult$freq,
      pdf_fd     = pdf_fd,
      logdensfd  = logdensfd,
      C          = C,
      denscdf    = denscdf,
      indcdf     = indcdf,
      Qvec       = Qvec,
      Hval       = Hval,
      DHval      = DHval,
      D2Hval     = D2Hval,
      active     = active,
      arclength  = arclength
    )
    
    parList[[icycle]] <- parListi
    
  }
  
  return(list(parList=parList))
  
} 