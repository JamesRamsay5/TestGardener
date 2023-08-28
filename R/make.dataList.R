
make.dataList <- function(U, key, optList, scrrng=NULL, titlestr=NULL,
                          nbin=nbinDefault(N), NumBasis=7, WfdPar=NULL,
                          jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
                          verbose=FALSE) {
  
  #' This function sets up the information required to analyze a set of data.
  #' The information is stored in the struct object dataStr.
  #' The information set up here is not meant to be modified by later code
  #' in the analysis.
  
  #  Last modified 15 May 2023 by Jim Ramsay
  
  N <- nrow(U)
  n <- ncol(U)
  
  # check U 
  
  if (min(U) < 1) stop(paste("Zero data values encountered in U.",
                             "Are they score values?"))
  Unumeric <- is.numeric(U)
  if (max(abs(floor(Unumeric)-Unumeric)) > 0) 
    stop("Non-integer values found in U.")
  
  #  check key
  
  if (!is.null(key)) {
    if (!is.numeric(key)) stop("Argument key is not numeric.")
    if (min(key) < 1) stop(paste("Zero data values encountered in key.",
                                 "Are they score values?"))
    if (max(abs(as.integer(key)-key)) > 0) 
      stop("Non-integer values found in key.")
   if (length(key) != n && length(key) > 0) 
      stop("length of key is neither n or 0") 
  }
  
  #  check optList
  
  if (!is.list(optList)) stop("Argument optList is not a list object.")
  if (length(optList) != 3) stop("Argument optList is not of length 3.")
  
  #  check scrrng
  
  if (!is.null(scrrng) & !is.numeric(scrrng))
    stop("Argument scrrng is neither NULL nor numeric.")
  if (is.numeric(scrrng) & length(scrrng) != 2) 
    stop("Argument scrrng is not of length 2.")
  
  ##  Set up noption using optList$optScr
  
  noption <- rep(0,n)
  for (item in 1:n) {
    noption[item] <- length(optList$optScr[[item]])  
  }
  
  #  construct logical vector grbg 
  
  grbg <- rep(FALSE, n)
  for (item in 1:n) {
    grbgind <- U[,item] > noption[item]
    if (sum(grbgind) > 0) {
      #  indices greater than noption[item] encountered
      #  add an option the these labels, set grbg value to TRUE
      #  change indices to updated noption[item]
      noption[item]   <- noption[item] + 1
      optScri <- optList$optScr[[item]]
      optScri <- c(optScri,0)
      optList$optScr[[item]] <- optScri
      grbg[item]      <- TRUE
      U[grbgind,item] <- noption[item]
    }
  }
  
  # compute dimension of ambient space
  
  Wdim  <- sum(noption) 
  
  ## compute sum scores for both examinees and items
  
  scrvec <- matrix(0,N,1)
  itmvec <- matrix(0,n,1)
  for (item in 1:n) {
    for (j in 1:N) {
      scorevec  <- optList$optScr[[item]]
      scoreij   <- scorevec[U[j,item]]
      scrvec[j] <- scrvec[j] + scoreij
      itmvec[item] <- itmvec[item] + scoreij
    }
  }
  
  scrmin  <- min(scrvec)
  scrmax  <- max(scrvec)
  if (is.null(scrrng)) scrrng <- c(scrmin,scrmax)
  nfine   <- 101
  scrfine <- seq(scrrng[1],scrrng[2],len=nfine)
    
  thetaQnt <- seq(0,100,len=2*nbin+1)
  
  ##  jitter sum scores if required
  
  if (jitterwrd) {
    scrjit <- scrvec + rnorm(N)*0.1
    scrjit[scrjit < scrmin] <- scrmin
    scrjit[scrjit > scrmax] <- scrmax
  } else {
    scrjit = scan("scrjit.txt",0)
  }
  
  ##  compute ranks for jittered sum scores
  
  scrrnk <- matrix(0,N,1)
  for (j in 1:N) scrrnk[j] <- sum(scrjit <= scrjit[j])
  percntrnk <- 100*scrrnk/N
  
  ##  Basis and bin setup for W function and theta estimation cycle
  
  #  number of basis functions.  If NULL, this is assigned according to size of N
  
  if (is.null(WfdPar)) {
    #  Default fdPar objects for representing functions
    if (round(NumBasis) == NumBasis) {
      #  NumBasis is an integer
      if (NumBasis < 2) 
        stop("There must be at least two spline basis functions.")
      if (NumBasis > 7) 
        warning("More than 7 basis functions may cause instability 
                in optimization.")
      Wnbasis <- NumBasis
      Wnorder <- min(Wnbasis, 5)
      Wbasis  <- fda::create.bspline.basis(c(0,100), Wnbasis, Wnorder) 
      WfdPar  <- fdPar(Wbasis)
    } else {
      stop("Number of basis functions is not an integer.")
    }
  } else {
    if (inherits(WfdPar, "basisfd") || inherits(WfdPar, "fd"))
      WfdPar <- fdPar(WfdPar)
    if (!inherits(WfdPar, "fdPar")) stop("WfdPar is not an fdPar object.")
  }
  
  ##  Wbinsmth.init computes an approximation to optimal Bmat
  
  WfdList <- Wbinsmth.init(percntrnk, nbin, WfdPar, grbg, optList, U) 
  
  ##  Construct dataList object to define data Listucture
  
  dataList <- list(U           = U, 
                   N           = N,
                   n           = n,
                   titlestr    = titlestr,
                   optList     = optList,
                   WfdList     = WfdList,
                   key         = key,
                   grbg        = grbg,
                   WfdPar      = WfdPar, 
                   noption     = noption, 
                   nbin        = nbin, 
                   scrrng      = scrrng, 
                   scrfine     = scrfine,
                   scrvec      = scrvec,
                   scrjit      = scrjit,
                   itmvec      = itmvec, 
                   percntrnk   = percntrnk, 
                   thetaQnt    = thetaQnt,
                   Wdim        = Wdim, 
                   PcntMarkers = PcntMarkers,
                   titlestr    = titlestr,
                   grbg        = grbg)
  
  return(dataList)
  
}

#  ---------------------------------------------------------------

Wbinsmth.init <- function(percntrnk, nbin, WfdPar, grbg, optList, U) {
  
  # Last modified 19 December 2022 by Jim Ramsay
  
  #  This version of Wbinsmth() uses direct least squares smoothing of the
  #  surprisal values at bin centers to generate dependent variables for
  #  a model for the vectorized K by M-1 parameter matrix Bmat.
  
  nitem <- ncol(U)
  chartList <- vector("list", nitem)
  indfine   <- seq(0,100, len=101)
  thetaQnt  <- seq(0,100, len=2*nbin+1)  
  bdry      <- thetaQnt[seq(1,2*nbin+1,by=2)]
  aves      <- thetaQnt[seq(2,2*nbin,  by=2)]  
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
    
    #  generate a map into M-vectors with zero row sums
    
    if (Mi == 2) {
      root2 <- sqrt(2)
      Zmati <- matrix(1/c(root2,-root2),2,1)
    } else {
      Zmati <- fda::zerobasis(Mi)
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
    Bmati  <- Result$coefficient
    Wfdi   <- fda::fd(Bmati, Wbasis)
    
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

#  ---------------------------------------------------------------

nbinDefault <- function(N) {
  if (N <= 500)              nbin <- floor(N/25)  
  if (N >  500 && N <= 2000) nbin <- floor(N/50)  
  if (N > 2000 && N <= 1e4)  nbin <- floor(N/100) 
  if (N >  1e4)              nbin <- 100 
  return(nbin)
}
  



    
    
