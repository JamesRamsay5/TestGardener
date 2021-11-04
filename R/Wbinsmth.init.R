Wbinsmth.init <- function(percntrnk, nbin, WfdPar, grbg, optList, U) {
  
  # Last modified 2 November 2021 by Jim Ramsay
  
  #  This version of Wbinsmth() uses direct least squares smoothing of the
  #  surprisal values at bin centers to generate dependent variables for
  #  a linear model for the vectorized K by M-1 parameter matrix Bmat.
  
  nitem <- ncol(U)
  chartList <- vector("list", nitem)
  indfine  <- seq(0,100, len=101)
  thetaQnt <- seq(0,100, len=2*nbin+1)  
  bdry     <- thetaQnt[seq(1,2*nbin+1,by=2)]
  aves     <- thetaQnt[seq(2,2*nbin,  by=2)]  
  freq <- matrix(0,nbin,1)
  freq[1] <- sum(percntrnk < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < percntrnk & percntrnk <= bdry[k])
  }
  meanfreq <- mean(freq)
  WfdList  <- vector("list",nitem)
  Wfd      <- WfdPar$fd
  Wbasis   <- Wfd$basis
  Wnbasis  <- Wbasis$nbasis
  for (item in 1:nitem) {
    Mi    <- length(optList$optScr[[item]])
    logMi <- log(Mi)
    Uveci <- as.numeric(U[,item])
    Pbin  <- matrix(0,nbin,Mi)  #  probabilities
    Wbin  <- matrix(0,nbin,Mi)  #  transformation of probability
    for (k in 1:nbin) {
      #  index of percntrnk values within this bin
      indk   <- percntrnk >= bdry[k] & percntrnk <= bdry[k+1]
      if (sum(indk) > 0) {
        Uvecik <- Uveci[indk]
        nk     <- sum(indk)
        for (m in 1:Mi) {
          Pbin[k,m] <- sum(Uvecik == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
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
      if (m != grbg[item]) {
        #  normal non-garbage choice, change NA values to SurprisalMax
        Wmis.na <- is.na(Pbin[,m])
        Wbin[Wmis.na,m] <- SurprisalMax
      } else {
        #  garbage choices: compute sparse numeric values into 
        #  linear approximationsand NA values to SurprisalMax
        Wmis.na <- is.na(Pbin[,m])
        indm <- (1:nbin)[!Wmis.na]
        indmlen <- length(indm)
        if (indmlen > 3) {
          WY <- Wbin[indm,m];
          WX <- cbind(rep(1,indmlen), aves[indm])
          BX <- lsfit(aves[indm], WY)$coefficients
          Wbin[indm,m] <- WX %*% BX
          Wbin[Wmis.na,m] <- SurprisalMax
        } else {
          Wbin[Wmis.na,m] <- SurprisalMax
        }
      }
    }
    
    #  generate a map into M-vectors with zero row sums
    if (Mi == 2) {
      root2 <- sqrt(2)
      Zmati <- matrix(1/c(root2,-root2),2,1)
    } else {
      Zmati <- zerobasis(Mi)
    }
    
    #  apply conventional smoothing of surprisal values
    Sfdi     <- fda::smooth.basis(aves, Wbin, WfdPar)$fd
    #  compute spline basis functions at bin centres
    Phimati  <- fda::eval.basis(aves, WfdPar$fd$basis)
    #  evaluate smooth at bin centres
    Smathati <- fda::eval.fd(aves, Sfdi)
    #  map this into zero-row-sum space
    Smatctri <- Smathati %*% Zmati
    #  regress the centred data on the negative of basis values
    Result <- lsfit(-Phimati, Smatctri, intercept=FALSE)
    Bmati  <- Result$coeff
    Wfdi   <- fd(Bmati, Wbasis)
    
    #  store objects in WListi
    
    WListi <- list(
      Wfd        = Wfdi,       #  functional data object for (options
      M          = Mi,         #  the number of options
      Pbin       = Pbin,       # proportions at each bin
      Wbin       = Wbin,       # negative surprisals at each bin
      Pmatfine   = NULL,   
      Wmatfine   = NULL,   
      DWmatfine  = NULL,  
      D2Wmatfine = NULL  
    )
    WfdList[[item]] <- WListi
  }
  
  return(WfdList)
  
}
