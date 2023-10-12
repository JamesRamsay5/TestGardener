Wpca <- function(WfdList, Wdim=NULL, nharm=2, rotate=TRUE) {
  
  #  Last modified 10 October 2023 by Jim Ramsay
  
  #  set up the dimension of the over-space containing the test info curve
  
  if (is.null(Wdim)) {
    Wdim <- 0
    for (i in 1:length(WfdList)) {
      WfdListi <- WfdList[[i]]
      Wdim <- Wdim + WfdListi$M
    }
  }
  
  #  compute grid indfine the probability and surprisal values,
  #  and the first derivative of the surprisal values for the total 
  #  of the curves
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  Pmat_full  <- matrix(0,nfine, Wdim)
  Wmat_full  <- matrix(0,nfine, Wdim)
  DWmat_full <- matrix(0,nfine, Wdim)
  
  n <- length(WfdList)
  m2 <- 0
  for (item in 1:n) {
    WListi <- WfdList[[item]]
    Wfdi   <- WListi$Wfd
    Mi     <- WListi$M
    m1 <- m2 +  1
    m2 <- m2 + Mi
    Wmat_full[,m1:m2]  <- eval.surp(indfine, Wfdi)
    DWmat_full[,m1:m2] <- eval.surp(indfine, Wfdi,1)
    Pmat_full[,m1:m2]  <- Mi^(Wmat_full[,m1:m2])
  }
  
  #  set up basis for (smoothing over [0,100])
  
  Wnbasis <- 7
  Wnorder <- 4
  Wbasis  <- fda::create.bspline.basis(c(0,100), Wnbasis, Wnorder) 
  WfdPar  <- fda::fdPar(fd(matrix(0,Wnbasis,1),Wbasis))
  
  #  smooth all Wdim surprisal curves to get best fitting fd curves
  
  Sfd <- fda::smooth.basis(indfine, Wmat_full, WfdPar)$fd
  
  #  functional PCA of the fd versions of surprisal curves
  #  the output is nharm principal component functions
  
  pcaList <- fda::pca.fd(Sfd, nharm, WfdPar, FALSE)
  
  #  set up the unrotated harmonic functional data object
  
  harmfd <- pcaList$harmonics
  
  #  compute variance proportions for unrotated solution
  
  varprop  <- pcaList$varprop
  if (!rotate) {
    varmxList    <- pcaList
    harmvarmxfd  <- harmfd
    varpropvarmx <- varprop
  } else {
    varmxList    <- varmx.pca.fd(pcaList)
    harmvarmxfd  <- varmxList$harmonics
    varpropvarmx <- varmxList$varprop
  }
  
  return(list(harmvarmxfd=harmvarmxfd, varpropvarmx=varpropvarmx))
  
}