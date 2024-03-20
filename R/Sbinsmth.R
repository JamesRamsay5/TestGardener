Sbinsmth <- function(index, dataList,  
                     indexQnt=seq(0,100,len=2*nbin+1), wtvec=matrix(1,n,1),
                     iterlim=20, conv=1e-4, dbglev=0) {
  # Compute smooth surprisal curves for binned data
  
  #  Arguments:
  #  index    ... vector of length N of score index values
  #  dataList ... List vector containing specifications of objects required
  #               to make an analysis and computed from input data
  #  indexQnt ... Sequence of boundary - binctr pairs followed by final boundary
  #  wtvec    ... Possible weights for items, but ordinarily all ones
  #  iterlim  ... Maximum number of iterations for optimising curves
  #  conv     ... Convergence criterion for iterations
  #  dbglev   ... Amoount of printed history of optimizations. None if zero
  
  # Last modified 15 January 2024 by Jim Ramsay

  #  -----------------------------------------------------------------------------
  #  Step 1.       Set up  objects required for subsequent steps
  #  -----------------------------------------------------------------------------

  #  check dataList
  
  if (!inherits(dataList,'list'))
    stop("Argument dataList is not of class list.")
  
  #  objects from dataList
  
  SfdList <- dataList$SfdList
  nbin    <- dataList$nbin
  Sbasis  <- dataList$Sbasis
  chcemat <- dataList$chcemat
  noption <- dataList$noption
  grbgvec <- dataList$grbgvec
  key     <- dataList$key

  # number of items plus at most one extra for guessing choices
  
  n     <- ncol(chcemat)
  nitem <- length(SfdList)  
  
  #  check objects from dataList
  
  if (is.null(chcemat)) stop("chcemat is null.") 
  if (is.null(noption)) stop("noption is null.")
  if (is.null(nbin))    stop("nbin is null.")
  
  #  -----------------------------------------------------------------------------
  #  Step 2. Bin the locations in index into bins defined by the
  #          bin edge and boundary vector indexQnt
  #  -----------------------------------------------------------------------------

  #  bin centers
  
  binctr <- indexQnt[seq(2,2*nbin,by=2)]
  
  #  bin boundaries, set at the corresponding quantiles in indexQnt.
  
  bdry <- indexQnt[seq(1,2*nbin+1,by=2)]
  bdry[1]      <-   0
  bdry[nbin+1] <- 100
  
  indfine <- seq(0,100,len=101)
  
  #  -----------------------------------------------------------------------------
  #  Step 3A.  Loop through the items to define the negative-surprisal S-curves
  #           for each question.
  #  -----------------------------------------------------------------------------
  
  for (item in 1:n) {

        SListi <- SfdList[[item]]
    Mi     <- SListi$M
    logMi  <- log(Mi)
    Sfdi   <- SListi$Sfd
    Zmati  <- SListi$Zmat
    
    #  set up bin probabilities and surprisals
    
    # print("entering binSP on line 80")
    result   <- binSP(item, index, dataList, bdry) 
    Pbin     <- result$Pbin
    Sbin     <- result$Sbin
    freq     <- result$freq
    meanfreq <- mean(freq)
    
    #  Set up SurprisalMax to replace NA's

    # print("entering surp.max on line 88")
    result <- surp.max(item, Mi, logMi, nbin, binctr, Pbin, Sbin, meanfreq, grbgvec)
    Pbin   <- result$Pbin
    Sbin   <- result$Sbin
    SurprisalMax <- result$SurprisalMax
    for (m in 1:Mi) {
      indNA <- is.na(Sbin[,m])
      Sbin[indNA,m] <- SurprisalMax
    }
    
    #  --------------------------------------------------------------------
    #  apply surprisal smoothing
    #  --------------------------------------------------------------------
    
    Bmat0   <- Sfdi$coefs
    Sbasis  <- Sfdi$basis
    #  suprisal smooth
    # if (outputwrd) dbglev <- 1 else dbglev <- 0
    #  optimize fit of smooth surprisal curves
    # print("entering smooth.surp on line 106")
    result  <- smooth.surp(binctr, Sbin, Bmat0, Sfdi, Zmati)
    #  retrieve results
    Sfdi    <- result$Sfd
    Bmati   <- result$Bmat
    # if (outputwrd) print(round(Zmati %*% t(Bmati),1))
    SSE     <- result$SSE
    hmat    <- result$hmat
    DvecSmatDvecB <- result$DvecSmatDvecB
    
    #  --------------------------------------------------------------------
    #  Step 4  Compute S and P values for bin point and each mesh point
    #  --------------------------------------------------------------------
    
    indfine   <- seq(0,100,len=101)
    # print("calling eval.surp line 120")
    Smatfine  <- eval.surp(indfine, Sfdi, Zmati)
    # print("calling eval.surp line 122")
    DSmatfine <- eval.surp(indfine, Sfdi, Zmati, 1)
    if (Sbasis$nbasis > 2) {
      # print("calling eval.surp line 125")
      D2Smatfine <- eval.surp(indfine, Sfdi, Zmati, 2)
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
                                     as.numeric(binctr))
    }
    #  compute derivatives for all options (except for diagonal)
    DPbinDS <- logMi*Pbinfit
    PStdErr <- sqrt(DPbinDS^2*SErrVar)
    SStdErr <- sqrt(SErrVar)  

    #  --------------------------------------------------------------------
    #  Step 7a assemble list vector SList that contains results
    #  of step 4 for this single item.  SfdList[[item]] <- SListi
    #  --------------------------------------------------------------------
    
    SListi  <- list(
      M          = Mi,         # the number of options
      Sfd        = Sfdi,       # functional data object for surprisal smooth
      Zmat       = Zmati,
      Pbin       = Pbin,       # proportions at each bin
      Sbin       = Sbin,       # negative surprisals at each bin
      indfine    = indfine,    # 101 equally spaced plotting points
      Pmatfine   = Pmatfine,   # Probabilities over fine mesh
      Smatfine   = Smatfine,   # S functions over fine mesh
      DSmatfine  = DSmatfine,  # 1st derivative of S functions over fine mesh
      D2Smatfine = D2Smatfine, # 2nd derivative of S functions over fine mesh
      PStdErr    = PStdErr,    # Std error for probabilities over fine mesh
      SStdErr    = SStdErr,    # Std error for surprisals    over fine mesh
      infoSurp   = infoSurp    # arc length of item information curve
    )

    SfdList[[item]] = SListi
    
  }
  
  # print("After loop through items")
  # print(class(SfdList))
  
  #  -----------------------------------------------------------------------------
  #  Step 3B.  Define the negative-surprisal S-curves for the extra question.
  #  -----------------------------------------------------------------------------
  
  if (nitem == n+1) {
    
    #  --------------------------------------------------------------------
    #  Step 3b assemble list vector SList that contains results
    #  of step 4 for this single item.  SfdList[[nitem]] <- SListi
    #  --------------------------------------------------------------------
    
    SListi <- SfdList[[nitem]]
    Mi    <- 3
    logMi <- 1.098612
    Sfdi  <- SListi$Sfd
    Zmati <- t(matrix(c( 0.7071068,  0.4082483, 0.0, 
                         -0.8164966, -0.7071068, 0.4082483),2,3))
    
    Pbin  <- matrix(0,nbin,Mi)  #  probabilities
    Sbin  <- matrix(0,nbin,Mi)  #  transformation of probability
    
    for (k in 1:nbin) {
      if (k == 1) {
        indk <- index <= bdry[2] 
      } else {
        indk <- index > bdry[k] & index <= bdry[k+1]
      }
      for (item in 1:n) {
        keyi <- key[item]
        chceveci <- as.numeric(chcemat[,item])
        indi1 <- chceveci != keyi & chceveci == 1
        Pbin[k,1] <- Pbin[k,1] + sum(indk & indi1)
        indi3 <- chceveci != keyi & chceveci == 3
        Pbin[k,2] <- Pbin[k,2] + sum(indk & indi3)
        # probability conditional on not choosing right answer
        # as well as not 1 or 3
        # indix <- chceveci != keyi & chceveci != 1 & chceveci != 3
        # probability conditional on not choosing 1 or 3
        indix <- chceveci != 1 & chceveci != 3
        Pbin[k,3] <- Pbin[k,3] + sum(indk & indix)
      }
      Pbinsumk <- sum(Pbin[k,])
      Pbin[k,] <- Pbin[k,]/Pbinsumk
      Sbin[k,] <- -log(Pbin[k,])/logMi
    }
    for (m in 1:Mi) {
      indInf <- is.infinite(Sbin[,m])
      Sbin[indInf,m] <- SurprisalMax
    }
    
    Bmat0   <- Sfdi$coefs
    Sbasis  <- Sfdi$basis
    #  suprisal smooth
    # if (outputwrd) dbglev <- 1 else dbglev <- 0
    #  optimize fit of smooth surprisal curves
    result  <- smooth.surp(binctr, Sbin, Bmat0, Sfdi, Zmati)
    #  retrieve results
    Sfdi    <- result$Sfd
    Bmati   <- result$Bmat
    # if (outputwrd) print(round(Zmati %*% t(Bmati),1))
    SSE     <- result$SSE
    hmat    <- result$hmat
    DvecSmatDvecB <- result$DvecSmatDvecB
    
    indfine   <- seq(0,100,len=101)
    # print("calling eval.surp line 262")
    Smatfine  <- eval.surp(indfine, Sfdi, Zmati)
    # print("calling eval.surp line 264")
    DSmatfine <- eval.surp(indfine, Sfdi, Zmati, 1)
    if (Sbasis$nbasis > 2) {
      # print("calling eval.surp line 267")
      D2Smatfine <- eval.surp(indfine, Sfdi, Zmati, 2)
    } else {
      D2Smatfine <- NULL
    }
    Pmatfine   <- Mi^(-Smatfine)
    
    infofine <- pracma::cumtrapz(indfine,sqrt(apply(DSmatfine^2,1,sum)))
    infoSurp <- infofine[101]
    
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
                                     as.numeric(binctr))
    }
    #  compute derivatives for all options (except for diagonal)
    DPbinDS <- logMi*Pbinfit
    PStdErr <- sqrt(DPbinDS^2*SErrVar)
    SStdErr <- sqrt(SErrVar)  
    
    #  --------------------------------------------------------------------
    #  Step 7 assemble list vector SList that contains results
    #  of step 4 for this single item.  SfdList[[ntem]] <- SListi
    #  --------------------------------------------------------------------
    
    SListi  <- list(
      M          = Mi,         # the number of options
      Sfd        = Sfdi,       # functional data object for surprisal smooth
      Zmat       = Zmati,
      Pbin       = Pbin,       # proportions at each bin
      Sbin       = Sbin,       # negative surprisals at each bin
      indfine    = indfine,    # 101 equally spaced plotting points
      Pmatfine   = Pmatfine,   # Probabilities over fine mesh
      Smatfine   = Smatfine,   # S functions over fine mesh
      DSmatfine  = DSmatfine,  # 1st derivative of S functions over fine mesh
      D2Smatfine = D2Smatfine, # 2nd derivative of S functions over fine mesh
      PStdErr    = PStdErr,    # Probabilities over fine mesh
      SStdErr    = SStdErr,    # S functions over fine mesh
      infoSurp   = infoSurp    # arc length of item information curve
    )
    
    SfdList[[nitem]] = SListi
    
    # print("After extra item")
    # print(class(SfdList))
    
  }
  
  return(list(SfdList=SfdList, binctr=binctr, bdry=bdry, freq=freq))
}

#. ----------------------------------------------------------------------------

binSP <- function(item, index, dataList, bdry) {
  
  #. Arguments
  
  #  index ... Vector of length N of score index values within [0,100]
  #. dataList ... List vector of length n or n+1 containing lists of information
  #               provided by make_dataList before analysis.
  #  bdry     ... Vector of length nbin + 1 of bin boundaries in [0,100]
  
  #. last modified 4 January 2024 by Jim
  
  #.  ---------------------------------------------------------------------------
  #             Bin the data and convert proportions to surpisal
  #.  ---------------------------------------------------------------------------
  
  nbin    <- dataList$nbin
  chcemat <- dataList$chcemat
  n       <- ncol(chcemat)
  SfdList <- dataList$SfdList
  
  #  compute frequencies for each bin
  
  freq <- matrix(0,nbin,1)
  freq[1] <- sum(index < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < index & index <= bdry[k])
  }
  
  meanfreq <- mean(freq)
  
  #  set some variable values for this item
  SListi   <- SfdList[[item]]
  Mi       <- SListi$M
  logMi <- log(Mi)
  #  extract active cases for (this item and corresponding index value
  chceveci <- as.numeric(chcemat[,item])
  #  bin frequencies (bin number nbin + 1 is the upper boundary)
  #  set up matrices to hold bin P and S values for each item and option
  Pbin <- matrix(0,nbin,Mi)  #  probabilities
  Sbin <- matrix(0,nbin,Mi)  #  transformation of probability
  #  --------------------------------------------------------------------
  #  Loop through the bins to compute binned P values and their
  #  transformation(s) to binned S values
  #  --------------------------------------------------------------------
  for (k in 1:nbin) {
    #  index of index values within this bin
    if (k == 1) {
      indk <- index <= bdry[2]
    } else {
      indk <- index > bdry[k] & index <= bdry[k+1]
    }
    if (sum(indk) > 0) {
      #  ------------------------------------------------------------
      #   Compute P values
      #  ------------------------------------------------------------
      nk     <- sum(indk)
      for (m in 1:Mi) {
        #  item choice probabilities
        Pbin[k,m] <- sum(chceveci[indk] == m)/nk
        if (Pbin[k,m] == 0) Pbin[k,m] <- NA
      }
      #  ------------------------------------------------------------
      #  convert the values of Pbin to the values of Sbin
      #  ------------------------------------------------------------
      Sbin[k,] <- -log(Pbin[k,])/logMi
    } else {
      Pbin[k,] <- NA
      Sbin[k,] <- NA
    }
  } # end of bin loop
  
  return(list(Pbin=Pbin, Sbin=Sbin, freq=freq))
}

#. ----------------------------------------------------------------------------

surp.max <- function(item, Mi, logMi, nbin, binctr, Pbin, Sbin, 
                     meanfreq, grbgvec) {
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
        SX <- cbind(rep(1,indmlen), binctr[indm])
        BX <- lsfit(binctr[indm], SY)$coefficients
        Sbin[indm,m]    <- SX %*% BX
        Sbin[nonindm,m] <- SurprisalMax
      } else {
        Sbin[nonindm,m] <- SurprisalMax
      }
    }
  }
  return(list(Pbin=Pbin, Sbin=Sbin, SurprisalMax=SurprisalMax))
}

