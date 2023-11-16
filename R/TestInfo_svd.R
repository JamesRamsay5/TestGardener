TestInfo_svd <- function(scrfine, SfdList, itemindex=1:n, nharm=2) {
  
  #  Construct a singular value approximation of the test info  
  #  of dimensions.
  
  #  Last modified 31 October 2023 by Jim Ramsay
  
  n = length(SfdList)
  
  #  set up the dimension of the over-space containing the curve
  
  dim <- 0
  for (i in itemindex) dim <- dim + SfdList[[i]]$M
  
  #  compute the full Test Info curve along a fine mesh
  
  mat_full <- matrix(0,100, dim)
  m2   <- 0
  for (i in itemindex) {
    SfdListi <- SfdList[[i]]
    Mi       <- SfdListi$M
    m1       <- m2 +  1
    m2       <- m2 + Mi
    mat_full[,m1:m2] <- SfdListi$matfine
  }
  
  #  singular value decomposition of curve values
  
  svd_result <- svd(mat_full, dim, dim)
  ind  <- 1:nharm
  Uhat <- svd_result$u[,ind]
  Vhat <- svd_result$v[,ind]
  Dhat <- svd_result$d[ind]
  
  # construct appoximation of test info curve of dimension nharm
  
  curvehat <- Uhat %*% diag(Dhat) %*% t(Vhat)
  
  #. proportions of shape variation accounted for by dimensions
  
  propvar = 100*Dhat^2/sum(svd_result$d^2)
  
  #  pull curvehat into the surprisal surface and output as
  #  surprisal functional data object
  
  Sbasis = create.bspline.basis(c(min(scrfine),max(scrfine)), 7)
  STfd   = smooth.surp(scrfine, curvehat, matrix(0,7,2), Sbasis)
  
  return(curvehat)
}
  