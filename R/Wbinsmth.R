Wbinsmth <- function(theta, dataList, WfdList=dataList$WfdList, 
                     thetaQnt=seq(0,100, len=2*nbin+1), wtvec=matrix(1,n,1),
                     iterlim=20, conv=1e-4, dbglev=0) {
  
  # Last modified 24 April 2023 by Jim Ramsay

  #  -----------------------------------------------------------------------------
  #  Step 1.       Set up  objects required for subsequent steps
  #  -----------------------------------------------------------------------------

  #  check dataList
  
  if (!inherits(dataList,'list'))
    stop("Argument dataList is not of class list.")
  
  #  objects from dataList
  
  n       <- length(WfdList)
  indfine <- seq(0,100,len=101)
  nbin    <- dataList$nbin
  nitem   <- length(dataList$optList$optScr)
  WfdPar  <- dataList$WfdPar
  U       <- dataList$U
  noption <- dataList$noption
  grbg    <- dataList$grbg

  #  check objects from dataList
  
  if (is.null(U))       stop("U is null.") 
  if (is.null(indfine)) stop("indfine is null.")
  if (is.null(noption)) stop("noption is null.")
  if (is.null(nbin))    stop("nbin is null.")
  
  #  -----------------------------------------------------------------------------
  #  Step 2. Bin the locations in theta into bins defined by the
  #          bin edge and boundary vector thetaQnt
  #  -----------------------------------------------------------------------------

  #  bin boundaries, set at the corresponding quantiles in thetaQnt.
  
  bdry <- thetaQnt[seq(1,2*nbin+1,by=2)]
  bdry[1]      <-   0
  bdry[nbin+1] <- 100

  #  bin centers
  
  aves <- rep(0,nbin)
  for (k in 1:nbin) {
    aves[k] <- (bdry[k]+bdry[k+1])/2
  }

  #  compute frequencies for each bin
  
  freq <- matrix(0,nbin,1)
  freq[1] <- sum(theta < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < theta & theta <= bdry[k])
  }
  
  meanfreq <- mean(freq)
  
  #  set up objects for required defining WfdPar for each item
  
  Wbasis    <- WfdPar$fd$basis
  Wnbasis   <- Wbasis$nbasis
  WLfd      <- WfdPar$Lfd
  Wlambda   <- WfdPar$lambda
  Westimate <- WfdPar$estimate
  Wpenmat   <- WfdPar$penmat
  
  #  -----------------------------------------------------------------------------
  #  Step 3.  Loop through the items to define the negative-surprisal W-curves
  #           for each question.
  #  -----------------------------------------------------------------------------
  
  for (item in 1:nitem) {
    if (dbglev > 0) {
      print(paste("Item",item))
    }
    #  set some variable values for this item
    Mi    <- noption[item]
    logMi <- log(Mi)
    #  extract active cases for (this item and corresponding theta value
    Uveci     <- as.numeric(U[,item])
    #  bin frequencies (bin number nbin + 1 is the upper boundary)
    #  set up matrices to hold bin P and W values for each item and option
    Pbin <- matrix(0,nbin,Mi)  #  probabilities
    Wbin <- matrix(0,nbin,Mi)  #  transformation of probability
    #  --------------------------------------------------------------------
    #  Step 3.1  loop through the bins to compute binned P values and their
    #  transformation(s) to binned W values
    #  --------------------------------------------------------------------
    for (k in 1:nbin) {
      #  index of theta values within this bin
      indk   <- theta >= bdry[k] & theta <= bdry[k+1]
      if (sum(indk) > 0) {
        #  ------------------------------------------------------------
        #                 compute P values
        #  ------------------------------------------------------------
        Uvecik <- Uveci[indk]
        nk     <- sum(indk)
        for (m in 1:Mi) {
          Pbin[k,m] <- sum(Uvecik == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
        #  ------------------------------------------------------------
        #  convert the values of Pbin to the values of Wbin
        #  ------------------------------------------------------------
        Wbin[k,] <- -log(Pbin[k,])/logMi
      } else {
        Pbin[k,] <- NA
      }
    } # end of bin loop

    #  --------------------------------------------------------------------
    #  Step 3.2 Smooth the binned W values
    #  --------------------------------------------------------------------
    
    #  Set up SurprisalMax to replace NA's
    
    maxWbin <- 0
    for (m in 1:Mi) {
      Wmis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Wmis.na]
      if (length(indm) > 0) maxWbin <- max(c(maxWbin,max(Wbin[indm,m])))
    }
    SurprisalMax <- min(c(-log(1/(meanfreq*2))/logMi, maxWbin))
    
    #  process NA values in Wbin associated with zero probabilities
    
    for (m in 1:Mi) {
      Wmis.na <- is.na(Pbin[,m])
      if (!grbg[item] || (grbg[item] && m != Mi)) {
        Wbin[Wmis.na,m] <- SurprisalMax
      }  else {
        #  garbage choices: compute sparse numeric values into 
        #  linear approximations and NA values to SurprisalMax
        indm    <- (1:nbin)[!Wmis.na]
        indmlen <- length(indm)
        nonindm <- (1:nbin)[Wmis.na]
        if (indmlen > 3) {
          WY <- Wbin[indm,m];
          WX <- cbind(rep(1,indmlen), aves[indm])
          BX <- lsfit(aves[indm], WY)$coefficients
          Wbin[indm,m]    <- WX %*% BX
          Wbin[nonindm,m] <- SurprisalMax
        } else {
          Wbin[nonindm,m] <- SurprisalMax
        }
      }
    }
    
    #  apply surprisal smoothing
    
    WListi  <- WfdList[[item]]
    Wfdi    <- WListi$Wfd
    Bmat0   <- Wfdi$coefs
    result  <- smooth.surp(aves, Wbin, Bmat0, WfdPar)
    Wfdi    <- result$Wfd
    Bmati   <- result$Bmat
    SSE     <- result$SSE
    hmat    <- result$hmat
    DvecSmatDvecB <- result$DvecSmatDvecB
    
    #  --------------------------------------------------------------------
    #  Step 4  Compute W and P values for bin point and each mesh point
    #  --------------------------------------------------------------------
     
    Wmatfine   <- eval.surp(indfine, Wfdi)
    DWmatfine  <- eval.surp(indfine, Wfdi, 1)
    if (Wnbasis > 2) {
      D2Wmatfine <- eval.surp(indfine, Wfdi, 2)
    } else {
      D2Wmatfine <- NULL
    }
    Pmatfine   <- Mi^(-Wmatfine)
    
    #  ------------------------------------------------------------------------
    #  Step 5. Compute arc length values for equally spaced theta values.  
    #. Integration is by the trapezoidal rule.
    #  ------------------------------------------------------------------------
    
    arclengthvec <- pracma::cumtrapz(indfine,sqrt(apply(DWmatfine^2,1,sum)))
    arclength    <- arclengthvec[101]
    
    #  --------------------------------------------------------------------
    #  Step 6 Compute the standard errors of the W values for each 
    #  option.  These tend only to be used when plotting curves, and this 
    #  step can be made optional to improve speed.
    #  --------------------------------------------------------------------
      
    #  compute Wdf = trace(data2fitmat)
    #  SSE is squared errors summed across points and surprisal curves
    #  It is computed within smooth_surp by function surp_fit
    
    WResidVar <- as.numeric(SSE/(nbin*Mi))
    dataVar   <- DvecSmatDvecB %*% solve(hmat) %*% t(DvecSmatDvecB)
    WErrVar   <- WResidVar*matrix(diag(dataVar),nbin,Mi)
    #  compute the vector (binary case) or matrix (multi case) of
    #  the sampling variances for each bin (and each option)
    #  set up derivatives of P wrt W
    PStdErr <- matrix(0,nbin,Mi)
    WStdErr <- matrix(0,nbin,Mi)
    Pbinfit <- matrix(0,nbin,Mi)
    DPbinDW <- matrix(0,nbin,Mi)
    for (m in 1:Mi) {
      Pbinfit[,m] <- pracma::interp1(as.numeric(indfine), as.numeric(Pmatfine[,m]), 
                                     as.numeric(aves))
    }
    #  compute derivatives for all options (except for diagonal)
    DPbinDW <- logMi*Pbinfit
    PStdErr <- sqrt(DPbinDW^2*WErrVar)
    WStdErr <- sqrt(WErrVar)        
    
    #  --------------------------------------------------------------------
    #  Step 7 assemble list vector WList that contains results
    #  of step 4 for this single item.  WfdList[[item]] <- WListi
    #  --------------------------------------------------------------------
    
    WListi  <- list(
      Wfd        = Wfdi,       # functional data object 
      M          = Mi,         # the number of options
      Pbin       = Pbin,       # proportions at each bin
      Wbin       = Wbin,       # negative surprisals at each bin
      indfine    = indfine,    # 101 equally spaced plotting points
      Pmatfine   = Pmatfine,   # Probabilities over fine mesh
      Wmatfine   = Wmatfine,   # W functions over fine mesh
      DWmatfine  = DWmatfine,  # 1st derivative of W functions over fine mesh
      D2Wmatfine = D2Wmatfine, # 2nd derivative of W functions over fine mesh
      PStdErr    = PStdErr,    # Probabilities over fine mesh
      WStdErr    = WStdErr,    # W functions over fine mesh
      arclength  = arclength   # arc length of item information curve
    )

    WfdList[[item]] = WListi
    
  }

  return(list(WfdList=WfdList, aves=aves, bdry=bdry, freq=freq))

}

