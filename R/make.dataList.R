
make.dataList <- function(U, key, optList, grbg=rep(0,n), scrrng=NULL, titlestr=NULL,
                          nbin=nbinDefault(N), NumBasis=NULL, WfdPar=NULL,
                          jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
                          verbose=FALSE) {
  
#' This function sets up the information required to analyze a set of data.
#' The information is stored in the struct object dataStr.
#' The information set up here is not meant to be modified by later code
#' in the analysis.
  
#  Last modified 28 Sept 2021 by Jim Ramsay

N <- nrow(U)
n <- ncol(U)

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
    scoreij   <- optList$optScr[[i]][U[j,i]]
    if (length(scoreij) > 0) {
      if (is.na(scoreij))   print(paste("is.na score:",j))
      if (is.null(scoreij)) print(paste("is.null score:",j))
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
    scrjit <- scrvec
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
  #  FdPar object for representing functions
  #  The order of the B-splines is 5 because we need a 
  #  smooth first derivative.
  Wnorder <- 5  #  Order of the basis functions
  # Set up the basis object
  Wbasis  <- fda::create.bspline.basis(c(0,100), NumBasis, Wnorder) 
  Wlambda <- 1e4   #  smoothing parameter
  #  Compute the penalty matrix for the third derivative
  Wnderiv <- 3  
  Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
  #  Assemble this information into a fdPar object.
  WfdPar  <- fdPar(Wbasis, Wnderiv, Wlambda, TRUE, Wpenmat)
}

##  Set up initial WfdList. 

#  Wbinsmth.init computes an approximation to optimal Bmat

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
                 itmvec      = itmvec, 
                 percntrnk   = percntrnk, 
                 thetaQnt    = thetaQnt,
                 Wdim        = Wdim, 
                 PcntMarkers = PcntMarkers,
                 titlestr    = titlestr)

return(dataList)

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
  if (N >  500 && N <= 2000)  NumBasis <- round(-14.7 + 8*log10(N))  
  if (N > 2000 && N <= 1e4)   NumBasis <- round(-14.7 + 8*log10(N))  
  if (N >  1e4)               NumBasis <-  24                        
  return(NumBasis)
}



    
    
