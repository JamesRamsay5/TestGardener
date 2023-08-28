smooth.ICC <- function(x, item, theta, dataList, 
                       thetaQnt=seq(0,100, len=2*nbin+1), 
                       wtvec=matrix(1,n,1), iterlim=20, conv=1e-4, dbglev=0) {
  
  #  A version of Sbinsmth designed for a smoothing data for a single item
  #  defined by the contents of an ICC object.
  
  #  Arguments:
  #  x        ... an ICC object containing information about single item
  #  theta    ... A vector of length N of score index values within [0,100]
  #  dataList ... A list object containing initial set up information
  #  thetaQnt ... A vector of length 2*nbin+1
  #  wtvec    ... A vector of weights for values to be smoothed
  #  iterlim  ... Limit on iterations for smoothing function
  #  conv     ... Criterion for convergence for smoothing function
  #  dbglev   ... Output level during smoothing (0, 1, 2)
  
  # Last modified 7 August 2023 by Jim Ramsay

  ICC <- x
  
  #  -----------------------------------------------------------------------------
  #  Step 1.       Set up  objects required for subsequent steps
  #  -----------------------------------------------------------------------------

  #  check dataList
  
  if (!inherits(dataList,'list'))
    stop("Argument dataList is not of class list.")
  
  #  objects from dataList
  
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
  
    #  set some variable values for this item
    M    <- noption
    logM <- log(M)
    #  extract active cases for (this item and corresponding theta value
    Uvec     <- as.numeric(U[,item])
    #  bin frequencies (bin number nbin + 1 is the upper boundary)
    #  set up matrices to hold bin P and W values for each item and option
    Pbin <- matrix(0,nbin,M)  #  probabilities
    Sbin <- matrix(0,nbin,M)  #  transformation of probability
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
        Uveck <- Uvec[indk]
        nk     <- sum(indk)
        for (m in 1:M) {
          Pbin[k,m] <- sum(Uveck == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
        #  ------------------------------------------------------------
        #  convert the values of Pbin to the values of Sbin
        #  ------------------------------------------------------------
        Sbin[k,] <- -log(Pbin[k,])/logM
      } else {
        Pbin[k,] <- NA
      }
    } # end of bin loop

    #  --------------------------------------------------------------------
    #  Step 3.2 Smooth the binned W values
    #  --------------------------------------------------------------------
    
    #  Set up SurprisalMax to replace NA's
    
    maxSbin <- 0
    for (m in 1:M) {
      Smis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Smis.na]
      if (length(indm) > 0) maxSbin <- max(c(maxSbin,max(Sbin[indm,m])))
    }
    SurprisalMax <- min(c(-log(1/(meanfreq*2))/logM, maxSbin))
    
    #  process NA values in Sbin associated with zero probabilities
    
    for (m in 1:M) {
      Smis.na <- is.na(Pbin[,m])
      if (!grbg[item] || (grbg[item] && m != M)) {
        Sbin[Smis.na,m] <- SurprisalMax
      }  else {
        #  garbage choices: compute sparse numeric values into 
        #  linear approximations and NA values to SurprisalMax
        indm    <- (1:nbin)[!Smis.na]
        indmlen <- length(indm)
        nonindm <- (1:nbin)[Smis.na]
        if (indmlen > 3) {
          WY <- Sbin[indm,m];
          WX <- cbind(rep(1,indmlen), aves[indm])
          BX <- lsfit(aves[indm], WY)$coefficients
          Sbin[indm,m]    <- WX %*% BX
          Sbin[nonindm,m] <- SurprisalMax
        } else {
          Sbin[nonindm,m] <- SurprisalMax
        }
      }
    }
    
    #  apply surprisal smoothing
    
    Sfd    <- ICC$Wfd
    Bmat0  <- Sfd$coefs
    result <- smooth.surp(aves, Sbin, Bmat0, WfdPar)
    Sfd    <- result$Wfd
    Bmati  <- result$Bmat
    SSE    <- result$SSE
    hmat   <- result$hmat
    DvecSmatDvecB <- result$DvecSmatDvecB
    
    #  --------------------------------------------------------------------
    #  Step 4  Compute W and P values for bin point and each mesh point
    #  --------------------------------------------------------------------
     
    Sarray <- array(0, c(101,M,3))
    Sarray[,,1]   <- eval.surp(indfine, Sfd)
    Sarray[,,2]  <- eval.surp(indfine, Sfd, 1)
    if (Wnbasis > 2) {
      Sarray[,,3] <- eval.surp(indfine, Sfd, 2)
    } else {
      Sarray[,,3] <- NULL
    }
    Pmatfine   <- M^(-Sarray[,,1])
    
    #  --------------------------------------------------------------------
    #  Step 5  Compute the standard errors of the W values for each 
    #  option.  These tend only to be used when plotting curves, and this 
    #  step can be made optional to improve speed.
    #  --------------------------------------------------------------------
      
    #  compute Wdf = trace(data2fitmat)
    #  SSE is squared errors summed across points and surprisal curves
    #  It is computed within smooth_surp by function surp_fit
    
    WResidVar <- as.numeric(SSE/(nbin*M))
    dataVar   <- DvecSmatDvecB %*% solve(hmat) %*% t(DvecSmatDvecB)
    WErrVar   <- WResidVar*matrix(diag(dataVar),nbin,M)
    #  compute the vector (binary case) or matrix (multi case) of
    #  the sampling variances for each bin (and each option)
    #  set up derivatives of P wrt W
    PStdErr <- matrix(0,nbin,M)
    WStdErr <- matrix(0,nbin,M)
    Pbinfit <- matrix(0,nbin,M)
    DPbinDW <- matrix(0,nbin,M)
    for (m in 1:M) {
      Pbinfit[,m] <- pracma::interp1(as.numeric(indfine), as.numeric(Pmatfine[,m]), 
                                     as.numeric(aves))
    }
    #  compute derivatives for all options (except for diagonal)
    DPbinDW <- logM*Pbinfit
    PStdErr <- sqrt(DPbinDW^2*WErrVar)
    SStdErr <- sqrt(WErrVar)        
    
    #  --------------------------------------------------------------------
    #  Step 6 assemble list vector WList that contains results
    #  of step 4 for this single item.  WfdList[[item]] <- ICC
    #  --------------------------------------------------------------------
    
    ICC  <- list(
      Wfd        = Sfd,        # functional data object 
      M          = M,          # the number of options
      Pbin       = Pbin,       # proportions at each bin
      Sbin       = Sbin,       # negative surprisals at each bin
      indfine    = indfine,    # 101 equally spaced plotting points
      Sarray     = Sarray,     # Probabilities over fine mesh
      PStdErr    = PStdErr,    # Probabilities over fine mesh
      SStdErr    = SStdErr,    # S functions over fine mesh
      itemStr    = ICC$itemStr,  # Item label
      optStr     = ICC$optStr   # List vector of option labels
    )

  return(list(ICC=ICC, aves=aves, bdry=bdry, freq=freq))

}

