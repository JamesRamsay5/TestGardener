Sbinsmth <- function(index, dataList, SfdList=dataList$SfdList, 
                     indexQnt=seq(0,100, len=2*nbin+1), wtvec=matrix(1,n,1),
                     iterlim=20, conv=1e-4, dbglev=0) {
  
  # Last modified 3 November 2023 by Jim Ramsay

  #  -----------------------------------------------------------------------------
  #  Step 1.       Set up  objects required for subsequent steps
  #  -----------------------------------------------------------------------------

  #  check dataList
  
  if (!inherits(dataList,'list'))
    stop("Argument dataList is not of class list.")
  
  #  objects from dataList
  
  n       <- length(SfdList)
  indfine <- seq(0,100,len=101)
  nbin    <- dataList$nbin
  nitem   <- length(dataList$noption)
  SfdPar  <- dataList$SfdPar
  chcemat <- dataList$chcemat
  noption <- dataList$noption
  grbgvec <- dataList$grbgvec

  #  check objects from dataList
  
  if (is.null(chcemat))       stop("chcemat is null.") 
  if (is.null(indfine)) stop("indfine is null.")
  if (is.null(noption)) stop("noption is null.")
  if (is.null(nbin))    stop("nbin is null.")
  
  #  -----------------------------------------------------------------------------
  #  Step 2. Bin the locations in index into bins defined by the
  #          bin edge and boundary vector indexQnt
  #  -----------------------------------------------------------------------------

  #  bin boundaries, set at the corresponding quantiles in indexQnt.
  
  bdry <- indexQnt[seq(1,2*nbin+1,by=2)]
  bdry[1]      <-   0
  bdry[nbin+1] <- 100

  #  bin centers
  
  aves <- rep(0,nbin)
  for (k in 1:nbin) {
    aves[k] <- (bdry[k]+bdry[k+1])/2
  }

  #  compute frequencies for each bin
  
  freq <- matrix(0,nbin,1)
  freq[1] <- sum(index < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < index & index <= bdry[k])
  }
  
  meanfreq <- mean(freq)
  
  #  set up objects for required defining SfdPar for each item
  
  Sbasis    <- SfdPar$fd$basis
  Snbasis   <- Sbasis$nbasis
  SLfd      <- SfdPar$Lfd
  Slambda   <- SfdPar$lambda
  Sestimate <- SfdPar$estimate
  Spenmat   <- SfdPar$penmat
  
  #  -----------------------------------------------------------------------------
  #  Step 3.  Loop through the items to define the negative-surprisal S-curves
  #           for each question.
  #  -----------------------------------------------------------------------------
  
  for (item in 1:nitem) {
    if (dbglev > 0) {
      # print(paste("Item",item))
    }
    #  set some variable values for this item
    Mi    <- noption[item]
    logMi <- log(Mi)
    #  extract active cases for (this item and corresponding index value
    chceveci <- as.numeric(chcemat[,item])
    #  bin frequencies (bin number nbin + 1 is the upper boundary)
    #  set up matrices to hold bin P and S values for each item and option
    Pbin <- matrix(0,nbin,Mi)  #  probabilities
    Sbin <- matrix(0,nbin,Mi)  #  transformation of probability
    #  --------------------------------------------------------------------
    #  Step 3.1  loop through the bins to compute binned P values and their
    #  transformation(s) to binned S values
    #  --------------------------------------------------------------------
    for (k in 1:nbin) {
      #  index of index values within this bin
      indk   <- index >= bdry[k] & index <= bdry[k+1]
      if (sum(indk) > 0) {
        #  ------------------------------------------------------------
        #                 compute P values
        #  ------------------------------------------------------------
        nk     <- sum(indk)
        for (m in 1:Mi) {
          Pbin[k,m] <- sum(chceveci[indk] == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
        #  ------------------------------------------------------------
        #  convert the values of Pbin to the values of Sbin
        #  ------------------------------------------------------------
        Sbin[k,] <- -log(Pbin[k,])/logMi
      } else {
        Pbin[k,] <- NA
      }
    } # end of bin loop

    #  --------------------------------------------------------------------
    #  Step 3.2 Smooth the binned S values
    #  --------------------------------------------------------------------
    
    #  Set up SurprisalMax to replace NA's
    
    maxSbin <- 0
    for (m in 1:Mi) {
      Smis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Smis.na]
      if (length(indm) > 0) maxSbin <- max(c(maxSbin,max(Sbin[indm,m])))
    }
    SurprisalMax <- min(c(-log(1/(meanfreq*2))/logMi, maxSbin))
    
    #  process NA values in Sbin associated with zero probabilities
    
    for (m in 1:Mi) {
      Smis.na <- is.na(Pbin[,m])
      if (!grbgvec[item] || (grbgvec[item] && m != Mi)) {
        Sbin[Smis.na,m] <- SurprisalMax
      }  else {
        #  garbage choices: compute sparse numeric values into 
        #  linear approximations and NA values to SurprisalMax
        indm    <- (1:nbin)[!Smis.na]
        indmlen <- length(indm)
        nonindm <- (1:nbin)[Smis.na]
        if (indmlen > 3) {
          SY <- Sbin[indm,m];
          SX <- cbind(rep(1,indmlen), aves[indm])
          BX <- lsfit(aves[indm], SY)$coefficients
          Sbin[indm,m]    <- SX %*% BX
          Sbin[nonindm,m] <- SurprisalMax
        } else {
          Sbin[nonindm,m] <- SurprisalMax
        }
      }
    }
    
    #  apply surprisal smoothing
    
    SListi  <- SfdList[[item]]
    Sfdi    <- SListi$Sfd
    Bmat0   <- Sfdi$coefs
    result  <- smooth.surp(aves, Sbin, Bmat0, SfdPar)
    Sfdi    <- result$Sfd
    Bmati   <- result$Bmat
    SSE     <- result$SSE
    hmat    <- result$hmat
    DvecSmatDvecB <- result$DvecSmatDvecB
    
    #  --------------------------------------------------------------------
    #  Step 4  Compute S and P values for bin point and each mesh point
    #  --------------------------------------------------------------------
     
    Smatfine   <- eval.surp(indfine, Sfdi)
    DSmatfine  <- eval.surp(indfine, Sfdi, 1)
    if (Snbasis > 2) {
      D2Smatfine <- eval.surp(indfine, Sfdi, 2)
    } else {
      D2Smatfine <- NULL
    }
    Pmatfine   <- Mi^(-Smatfine)
    
    #  ------------------------------------------------------------------------
    #  Step 5. Compute arc length values for equally spaced index values.  
    #. Integration is by the trapezoidal rule.
    #  ------------------------------------------------------------------------
    
    infofine <- pracma::cumtrapz(indfine,sqrt(apply(DSmatfine^2,1,sum)))
    infoSurp <- infofine[101]
    
    #  --------------------------------------------------------------------
    #  Step 6 Compute the standard errors of the S values for each 
    #  option.  These tend only to be used when plotting curves, and this 
    #  step can be made optional to improve speed.
    #  --------------------------------------------------------------------
      
    #  compute Sdf = trace(data2fitmat)
    #  SSE is squared errors summed across points and surprisal curves
    #  It is computed within smooth_surp by function surp_fit
    
    SResidVar <- as.numeric(SSE/(nbin*Mi))
    dataVar   <- DvecSmatDvecB %*% solve(hmat) %*% t(DvecSmatDvecB)
    SErrVar   <- SResidVar*matrix(diag(dataVar),nbin,Mi)
    #  compute the vector (binary case) or matrix (multi case) of
    #  the sampling variances for each bin (and each option)
    #  set up derivatives of P wrt S
    PStdErr <- matrix(0,nbin,Mi)
    SStdErr <- matrix(0,nbin,Mi)
    Pbinfit <- matrix(0,nbin,Mi)
    DPbinDS <- matrix(0,nbin,Mi)
    for (m in 1:Mi) {
      Pbinfit[,m] <- pracma::interp1(as.numeric(indfine), as.numeric(Pmatfine[,m]), 
                                     as.numeric(aves))
    }
    #  compute derivatives for all options (except for diagonal)
    DPbinDS <- logMi*Pbinfit
    PStdErr <- sqrt(DPbinDS^2*SErrVar)
    SStdErr <- sqrt(SErrVar)        
    
    #  --------------------------------------------------------------------
    #  Step 7 assemble list vector SList that contains results
    #  of step 4 for this single item.  SfdList[[item]] <- SListi
    #  --------------------------------------------------------------------
    
    SListi  <- list(
      Sfd        = Sfdi,       # functional data object 
      M          = Mi,         # the number of options
      Pbin       = Pbin,       # proportions at each bin
      Sbin       = Sbin,       # negative surprisals at each bin
      indfine    = indfine,    # 101 equally spaced plotting points
      Pmatfine   = Pmatfine,   # Probabilities over fine mesh
      Smatfine   = Smatfine,   # S functions over fine mesh
      DSmatfine  = DSmatfine,  # 1st derivative of S functions over fine mesh
      D2Smatfine = D2Smatfine, # 2nd derivative of S functions over fine mesh
      PStdErr    = PStdErr,    # Probabilities over fine mesh
      SStdErr    = SStdErr,    # S functions over fine mesh
      infoSurp   = infoSurp   # arc length of item information curve
    )

    SfdList[[item]] = SListi
    
  }

  return(list(SfdList=SfdList, aves=aves, bdry=bdry, freq=freq))

}

