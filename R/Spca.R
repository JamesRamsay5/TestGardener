Spca <- function(SfdList, Sdim=NULL, nharm=2, rotate=TRUE) {
  
  #  Last modified 10 October 2023 by Jim Ramsay
  
  #  set up the dimension of the over-space containing the test info curve
  
  if (is.null(Sdim)) {
    Sdim <- 0
    for (i in 1:length(SfdList)) {
      SfdListi <- SfdList[[i]]
      Sdim <- Sdim + SfdListi$M
    }
  }
  
  #  compute grid indfine the probability and surprisal values,
  #  and the first derivative of the surprisal values for the total 
  #  of the curves
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  Pmat_full  <- matrix(0,nfine, Sdim)
  Smat_full  <- matrix(0,nfine, Sdim)
  DSmat_full <- matrix(0,nfine, Sdim)
  
  n <- length(SfdList)
  m2 <- 0
  for (item in 1:n) {
    SListi <- SfdList[[item]]
    Sfdi   <- SListi$Sfd
    Mi     <- SListi$M
    m1 <- m2 +  1
    m2 <- m2 + Mi
    Smat_full[,m1:m2]  <- eval.surp(indfine, Sfdi)
    DSmat_full[,m1:m2] <- eval.surp(indfine, Sfdi,1)
    Pmat_full[,m1:m2]  <- Mi^(Smat_full[,m1:m2])
  }
  
  #  set up basis for (smoothing over [0,100])
  
  Snbasis <- 7
  Snorder <- 4
  Sbasis  <- fda::create.bspline.basis(c(0,100), Snbasis, Snorder) 
  SfdPar  <- fda::fdPar(fd(matrix(0,Snbasis,1),Sbasis))
  
  #  smooth all Sdim surprisal curves to get best fitting fd curves
  
  Sfd <- fda::smooth.basis(indfine, Smat_full, SfdPar)$fd
  
  #  functional PCA of the fd versions of surprisal curves
  #  the output is nharm principal component functions
  
  pcaList <- fda::pca.fd(Sfd, nharm, SfdPar, FALSE)
  
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