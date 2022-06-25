
make.dataList <- function(U, key, optList, grbg=rep(0,n), scrrng=NULL, titlestr=NULL,
                          nbin=nbinDefault(N), NumBasis=NULL, WfdPar=NULL,
                          jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
                          quadwrd=FALSE, verbose=FALSE) {
  
  #' This function sets up the information required to analyze a set of data.
  #' The information is stored in the struct object dataStr.
  #' The information set up here is not meant to be modified by later code
  #' in the analysis.
  
  #  Last modified 13 June 2022 by Jim Ramsay
  
  N <- nrow(U)
  n <- ncol(U)
  
  # check U for having only positive integer values
  
  if (min(U) < 1) stop("Zero data values encountered in data.  Are they score values?")
  
  #  check lengths of key and garbg
  
  if (length(key) != n && length(key) > 0) {
    stop("length of key is neither n or 0") 
  } 
  
  if (length(grbg) != n && length(grbg) > 0) {
    stop("length of grbg is not n") 
  } 
  
  ##  Set up Wdim
  
  noption <- matrix(0,n,1)
  for (item in 1:n) {
    noption[item] <- length(optList$optScr[[item]])# JL 2021-02-18
  }
  Wdim   <- sum(noption) 
  
  ## compute sum scores for (both examinees and items
  
  scrvec <- matrix(0,N,1)
  itmvec <- matrix(0,n,1)
  for (i in 1:n) {
    for (j in 1:N) {
      scoreij <- optList$optScr[[i]][U[j,i]]
      if (!is.null(scoreij)) {
        if (is.na(scoreij)) print(paste("is.na score:",j))
      } else {
        print(paste("score of length 0:",j))
      }
      scrvec[j] <- scrvec[j] + scoreij
      itmvec[i] <- itmvec[i] + scoreij
    }
  }
  
  scrmin  <- min(scrvec)
  scrmax  <- max(scrvec)
  if (is.null(scrrng)) {
    scrrng <- c(scrmin,scrmax)
  }
  nfine   <- 101
  scrfine <- seq(scrrng[1],scrrng[2],len=nfine)
  denscdf <- seq(0,1,len=nfine)
  
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
  for (j in 1:N) {
    scrrnk[j] <- sum(scrjit <= scrjit[j])
  }
  
  percntrnk <- 100*scrrnk/N
  
  ##  Basis and bin setup for W function and theta estimation cycle
  
  #  number of basis functions.  If NULL, this is assigned according to size of N
  
  if (is.null(NumBasis)) {
    NumBasis=NumBasisDefault(N)
  }
  
  if (is.null(WfdPar)) {
    #  Defaualt fdPar objects for representing functions
    #  Two options are available: order 5 or order 3 (quadratic basis functions)
    if (quadwrd) {
      #  The order of the B-splines is 2 and the two basis functions are
      #  quad.  No smoothing is used.
      Wnorder <- 3  #  Order of the basis functions
      Wnbasis <- 3  #  Three basis functions
      # Set up the basis object
      Wbasis  <- fda::create.bspline.basis(c(0,100), Wnbasis, Wnorder) 
      WfdPar  <- fdPar(Wbasis)
    } else {
      #  The order of the B-splines is 5 because we need a 
      #  smooth first derivative.
      Wnorder <- 5  #  Order of the basis functions
      Wnbasis <- NumBasis  #  NumBasis basis functions
      # Set up the basis object
      Wbasis  <- fda::create.bspline.basis(c(0,100), Wnbasis, Wnorder) 
      Wlambda <- 1e4   #  smoothing parameter
      #  Compute the penalty matrix for the third derivative
      Wnderiv <- 3  
      Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
      #  Assemble this information into a fdPar object.
      WfdPar  <- fdPar(Wbasis, Wnderiv, Wlambda, TRUE, Wpenmat)
    } 
  }
  
  ##  Set up initial WfdList. 
  
  #  Wbinsmth.init computes an approximation to optimal Bmat
  
  # print("Entering Wbinsmth.init")
  # print(percntrnk[1:5])
  # print(nbin)
  # print(WfdPar$fd$basis$nbasis)
  
  WfdList <- Wbinsmth.init(percntrnk, nbin, WfdPar, grbg, optList, U) 
  
  ##  Construct dataList object to define data Listucture
  
  dataList <- list(U           = U, 
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
                   titlestr    = titlestr)
  
  return(dataList)
  
}

#  ---------------------------------------------------------------

Wbinsmth.init <- function(percntrnk, nbin, WfdPar, grbg, optList, U) {
  
  # Last modified 21 January 2022 by Jim Ramsay
  
  #  This version of Wbinsmth() uses direct least squares smoothing of the
  #  surprisal values at bin centers to generate dependent variables for
  #  a model for the vectorized K by M-1 parameter matrix Bmat.
  
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
    
    # if (item==1){
    #   print("Pbin 1:")
    #   print(Pbin)
    #   print("Wbin 1")
    #   print(Wbin)
    # }
    
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
    
    # if (item==1) print(paste("SurprisalMax",SurprisalMax))
          
    #  process NA values in Wbin associated with zero probabilities
    
    for (m in 1:Mi) {
      if (is.null(grbg[item]) || grbg[item] == 0) {
        Wmis.na         <- is.na(Pbin[,m])
        Wbin[Wmis.na,m] <- SurprisalMax
      }  else {
        if (m != grbg[item]) {
          #  normal non-garbage choice, change NA values to SurprisalMax
          Wmis.na         <- is.na(Pbin[,m])
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
    }
    
    # if (item==1) {
    #   print("Second Wbin 1")
    #   print(Wbin)
    # }
    
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
    Bmati  <- Result$coefficient
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

#  ---------------------------------------------------------------

nbinDefault <- function(N) {
  if (N <= 500)              nbin <- floor(N/25)  
  if (N >  500 && N <= 2000) nbin <- floor(N/50)  
  if (N > 2000 && N <= 1e4)  nbin <- floor(N/100) 
  if (N >  1e4)              nbin <- 100 
  return(nbin)
}
  
#  ---------------------------------------------------------------

NumBasisDefault <- function(N) {
  if (N <= 500)               NumBasis <- 7                          
  if (N >  500 && N <= 1e4)   NumBasis <- round(-14.7 + 8*log10(N))  
  if (N >  1e4)               NumBasis <-  24                        
  return(NumBasis)
}



    
    
