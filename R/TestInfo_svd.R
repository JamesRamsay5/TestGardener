TestInfo_svd <- function(scrfine, WfdList, itemindex=1:n, nharm=2) {
  
  #  Construct a singular value approximation of the test info  
  #  of dimensions.
  
  #  Last modified 5 April 2023 by Jim Ramsay
  
  n = length(WfdList)
  
  #  set up the dimension of the over-space containing the curve
  
  Wdim <- 0
  for (i in itemindex) Wdim <- Wdim + WfdList[[i]]$M
  
  #  compute the full Test Info curve along a fine mesh
  
  Wmat_full <- matrix(0,100, Wdim)
  m2   <- 0
  for (i in itemindex) {
    WfdListi <- WfdList[[i]]
    Mi       <- WfdListi$M
    m1       <- m2 +  1
    m2       <- m2 + Mi
    Wmat_full[,m1:m2] <- WfdListi$Wmatfine
  }
  
  #  singular value decomposition of curve values
  
  svd_result <- svd(Wmat_full, Wdim, Wdim)
  ind  <- 1:nharm
  Uhat <- svd_result$u[,ind]
  Vhat <- svd_result$v[,ind]
  Dhat <- svd_result$d[ind]
  
  # construct appoximation of test info curve of dimension nharm
  
  Wcurvehat <- Uhat %*% diag(Dhat) %*% t(Vhat)
  
  #. proportions of shape variation accounted for by dimensions
  
  propvar = 100*Dhat^2/sum(svd_result$d^2)
  
  #  pull Wcurvehat into the surprisal surface and output as
  #  surprisal functional data object
  
  Sbasis = create.bspline.basis(c(min(scrfine),max(scrfine)), 7)
  STfd   = smooth.surp(scrfine, Wcurvehat, matrix(0,7,2), Sbasis)
  
  
  return(Wcurvehat)
}
  