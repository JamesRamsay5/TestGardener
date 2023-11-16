Analyze <- function(index, indexQnt, dataList, ncycle = 10, itdisp = FALSE, 
                    verbose = FALSE) {
  
  # Last modified 2 November 2023 by Jim Ramsay

  #  set up list vector to contain all results for each cycle
  
  parmListvec <- vector("list",ncycle)  
  
  #  define the spline basis for representing the log density function
  
  logdensbasis <- fda::create.bspline.basis(c(0,100), 15)
  
  #  initialize surprisal curves and compute size of overspace
  
  SfdList <- dataList$SfdList
  
  n = length(SfdList)
  
  Sdim = 0
  for (i in 1:n) {
    SStri = SfdList[[i]]
    Sdim  = Sdim + SStri$M
  }
  
  #  set up data matrix
  
  chcemat <- dataList$chcemat
  
  #  main cycle loop
  
  for (icycle in 1:ncycle) {
    
    if (verbose) print(paste('----------  Cycle ',icycle,'-----------'))
    
    #  ----------------------------------------------------------
    #  Step 1:  Bin the data, and smooth the binned data
    #  ----------------------------------------------------------
    # # print("step 1")
    if (verbose) print("Optimize surprisal curves:")
    
    SfdResult <- Sbinsmth(index, dataList, SfdList, indexQnt)
    SfdList   <- SfdResult$SfdList
    binctr    <- SfdResult$aves
    bdry      <- SfdResult$bdry
    freq      <- SfdResult$freq
    
    # print("Step 1 finished")
    
    #  ----------------------------------------------------------
    #  Step 2:  compute mean value of objective function 
    #  ----------------------------------------------------------
    
    if (verbose) print("Compute mean examinee fits")
    
    Fvec  <- Ffun(index, SfdList, chcemat)
    meanF <- mean(Fvec)
    
    if (verbose) print(paste('Mean data fit = ', round(meanF,3)))
    # print("Step 2 finished")
    
    #  ----------------------------------------------------------
    #  Step 3:  Compute optimal score index values
    #  ----------------------------------------------------------
    # # print("step 3")
    
    if (verbose) print("Optimize examinee data fits")

    indexfunList <- index_fun(index, SfdList, chcemat, 20, 1e-3, itdisp=itdisp)
    index    <- indexfunList$index_out
    Fval     <- indexfunList$Fval
    DFval    <- indexfunList$DFval
    D2Fval   <- indexfunList$D2Fval
    active   <- indexfunList$active
    # print("Step e finished")
    
    #  ----------------------------------------------------------
    #  Step 4:  Estimate the score density for score index values
    #  ----------------------------------------------------------
    # # print("step 4")
    
    if (verbose) print("Compute score index density")
    
    indexdens <- index[0 < index & index < 100]
    index_distnList <- index_distn(indexdens, logdensbasis)
    pdf_fd    <- index_distnList$pdf_fd
    logdensfd <- index_distnList$logdensfd
    cdffine   <- index_distnList$cdffine
    C         <- index_distnList$C
    denscdf   <- index_distnList$denscdf
    indcdf    <- index_distnList$indcdf
    markers   <- dataList$PcntMarkers/100
    Qvec      <- pracma::interp1(as.numeric(denscdf), as.numeric(indcdf), markers)
    nbin      <- dataList$nbin
    bdry      <- seq(0,2*nbin,1)/(2*nbin)
    indexQnt  <- pracma::interp1(as.numeric(denscdf), as.numeric(indcdf), bdry)
    # print("Step 4 finished")
    
    #  ----------------------------------------------------------
    #  Step 5.  Compute arc length and its measures
    #  ----------------------------------------------------------
    
    DSfine = matrix(0,101,Sdim)
    m2 = 0
    for (i in 1:n) {
      SListi = SfdList[[i]]
      Mi     = SListi$M
      m1 = m2 + 1
      m2 = m2 + Mi
      DSfine[,m1:m2] = SListi$DSmatfine
    }
    indfine <- seq(0,100,len=101)
    infoSurp = max(pracma::cumtrapz(indfine, sqrt(apply(DSfine^2,1,sum))))
    
    if (verbose)  print(paste('infoSurp in bits = ',round(infoSurp,1)))
    # print("Step 5 finished")
    
    #  ----------------------------------------------------------
    #  Step 6:  Check for mis-identications of minimum index 
    #  ----------------------------------------------------------
    
    if (dataList$SfdPar$fd$basis$nbasis > 2) {
      Result <- index_search(SfdList, dataList$chcemat, index, Fval, DFval, D2Fval)
      index  <- Result$index
      Fval   <- Result$Fval
      DFval  <- Result$DFval
      D2Fval <- Result$D2Fval
    }
    # print("Step 6 finished")
    
    #  ----------------------------------------------------------
    #  Step 7:  set up ParameterCell arrays
    #  ----------------------------------------------------------

    parmListi <- list(
      index      = index,
      indexQnt   = indexQnt,
      SfdList    = SfdList,
      meanF      = meanF,
      binctr     = binctr,
      bdry       = bdry,
      freq       = freq,
      pdf_fd     = pdf_fd,
      logdensfd  = logdensfd,
      C          = C,
      denscdf    = denscdf,
      indcdf     = indcdf,
      Qvec       = Qvec,
      Fval       = Fval,
      DFval      = DFval,
      D2Fval     = D2Fval,
      active     = active,
      infoSurp   = infoSurp
    )
    
    parmListvec[[icycle]] <- parmListi
    # print("Step 7 finished")
    
  }
  
  return(parmListvec=parmListvec)
  
}