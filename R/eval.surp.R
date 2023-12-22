eval.surp <- function(evalarg, Sfdobj, Zmat, nderiv=0) {
  #  Evaluates a value of a surprisal coordinate functional data object. 
  #  Evaluates a value or a derivative of a surprisal coordinate functional  
  #  data object. 
  #  A positive functional data object h  is <- the form
  #           h(x) <- (exp Sfdobj)(x)
  #  Note that the first two arguments may be interchanged.
  #  
  #
  #  Arguments:
  #  EVALARG ... A vector of values at which all functions are to 
  #              evaluated.
  #  SFDOBJ  ... Functional data object.  It must define a single
  #              functional data observation.
  #  ZMAT    ... An M by M-1 matrix satisfying Z'Z = I and Z'1 = 0.
  #  NDERIV  ... The order of the derivative.  Must be 0, 1 or 2.
  #  Returns:  An array of function values corresponding to the 
  #              argument values in EVALARG
  
  #  Last modified 19 December 2023 by Jim Ramsay
  
  #  check arguments
  
  if (floor(nderiv) != nderiv) {
    stop('Third argument nderiv is not an integer.')
  }
  
  if (nderiv < 0 || nderiv > 2) {
    stop('Third argument nderiv is not 0, 1 or 2.')
  }
  
  #  Exchange the first two arguments if the first is an FD object
  #    and the second numeric
  
  if (is.numeric(Sfdobj) && fda::is.fd(evalarg)) {
    temp    <- Sfdobj
    Sfdobj  <- evalarg
    evalarg <- temp
  }
  
  #  Check the arguments
  
  if (!is.numeric(evalarg)) {
    stop('Argument EVALARG is not numeric.')
  }
  
  #  transpose EVALARG if necessary to make it a column vector
  
  evaldim <- dim(as.matrix(evalarg))
  if (evaldim[1] == 1 && evaldim[2] > 1) {  
    evalarg <- t(evalarg)  
  }
  
  #  check EVALARG 
  
  if (evaldim[1] > 1 && evaldim[2] > 1) 
    stop('Argument EVALARG is not a vector.')
    
  evalarg  <- as.vector(evalarg)
  
  #  check FDOBJ
  
  if (!fda::is.fd(Sfdobj)) 
    stop('Argument FD is not a functional data object.')
  
  #  check Sfdobj$basis
  
  basisfd  <- Sfdobj$basis
  rangeval <- basisfd$rangeval
  if (min(evalarg) < rangeval[1] || max(evalarg) > rangeval[2]) {
    stop("Values in argument valarg outside of basis boundaries.")
  }
  evalarg[evalarg < rangeval[1]-1e-10] <- NA
  evalarg[evalarg > rangeval[2]+1e-10] <- NA
  
  #  compute Zmat, a M by M-1 orthonormal matrix
  
  M     <- nrow(Zmat)
  Bmat  <- Sfdobj$coefs
  BZmat <- Bmat %*% t(Zmat)
  
  #  Set up coefficient array for FD
  
  Xmat     <- fda::eval.basis(evalarg, basisfd) %*% BZmat
  MXmat    <- M^Xmat
  sum0     <- rowSums(MXmat)
  SumMXmat <- matrix(rep(sum0,each=M), ncol=M, byrow=TRUE)
  
  #  Case where EVALARG is a vector of values to be used for all curves
  #  NB:  Resist the temptation of use /log(M) anywhere else.  
  #  This code sets up Smat as defined with log basis M, and no further
  #  definition of log basis is required.
  
  if (nderiv == 0) {
    Smat   <- -Xmat + log(SumMXmat)/log(M)
    return(Smat)
  }
  
  #  First derivative:
  
  if (nderiv == 1) {
    Pmat    <- MXmat/SumMXmat
    DXmat   <- fda::eval.basis(evalarg, basisfd, 1) %*% BZmat
    matSum  <- rowSums(Pmat*DXmat)
    Rmat    <- matrix(rep(matSum,each=M), ncol=M, byrow=TRUE)
    DSmat   <- -DXmat + Rmat
    return(DSmat)
  }
  
  #  Second derivative
  
  if (nderiv == 2) {
    Pmat    <- MXmat/SumMXmat
    DXmat   <- fda::eval.basis(evalarg, basisfd, 1) %*% BZmat
    D2Xmat  <- fda::eval.basis(evalarg, basisfd, 2) %*% BZmat
    matSum  <- rowSums(Pmat*DXmat)
    Rmat    <- matrix(rep(matSum,each=M), ncol=M, byrow=TRUE)
    matSum2 <- rowSums((Pmat*(D2Xmat + DXmat^2) - Rmat*DXmat))
    D2Smat  <- -D2Xmat + 
      matrix(rep(matSum2,each=M), ncol=M, byrow=TRUE)
    return(D2Smat)
  }
}

