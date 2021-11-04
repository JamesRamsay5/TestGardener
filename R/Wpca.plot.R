Wpca.plot <- function(arclength, WfdList, Wdim, nharm=2, rotate=TRUE, 
                       dodge = 1.003, titlestr=NULL) {
  
  #  Last modified 8 February 2021 by Jim Ramsay
  
  #  set up matrices of fine mesh values
  
  nfine   <- 101
  indfine <- seq(0,arclength,len=nfine)
  
  n <- length(WfdList)
  
  Pmat_full  <- matrix(0,nfine, Wdim)
  Wmat_full  <- matrix(0,nfine, Wdim)
  DWmat_full <- matrix(0,nfine, Wdim)
  
  m2 <- 0
  for (i in 1:n) {
    WListi <- WfdList[[i]]
    Wfdi   <- WListi$Wfd
    Mi     <- WListi$M
    m1 <- m2 +  1
    m2 <- m2 + Mi
    Wmat_full[,m1:m2]  <- eval.surp(indfine,Wfdi)
    DWmat_full[,m1:m2] <- eval.surp(indfine,Wfdi,1)
    Pmat_full[,m1:m2]  <- Mi^(Wmat_full[,m1:m2])
  }
  
  #  set up basis for (smoothing over arclength
  
  Wnbasis <- 24
  Wnorder <-  5
  Wbasis  <- fda::create.bspline.basis(c(0,arclength), Wnbasis, Wnorder) 
  Wlambda <- 1e4   #  smoothing parameter
  Wnderiv <- 3  
  Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
  WfdPar  <- fda::fdPar(Wbasis, Wnderiv, Wlambda, TRUE, Wpenmat)
  
  #  prepare surprisal curves
  #  fine mesh of points along manifold of length arclength
  nfine <- dim(Wmat_full)[1]
  arclengthfine <- seq(0,arclength,len=nfine)
  #  smooth all Wdim curves to get surprisal curves
  Sfd <- fda::smooth.basis(arclengthfine, Wmat_full, WfdPar)$fd
  
  #  functional PCA of the surprisal curves
  pcaList <- fda::pca.fd(Sfd, nharm, WfdPar, FALSE)
  #  set up the unrotated harmonic functional data object
  harmfd <- pcaList$harmonics
  #  compute variance proportions for (unrotated solution
  varprop <- pcaList$varprop
  
  #  carry out a varimax rotation
  
  if (rotate) {
    varmxList    <- pcaList
    harmvarmxfd  <- harmfd
    varpropvarmx <- varprop
  } else {
    varmxList  <- varmx.pca.fd(pcaList)
    harmvarmxfd  <- varmxList$harmfd
    varpropvarmx <- varmxList$varprop
  }
  
  #  set up the varimax rotated harmonic functional data object
  #  display variance proportions for (varimax rotated solution
  
  #  plot the first two or three harmonics
  #  ggpot version
  if (nharm == 2 || nharm == 3) {
    pind   <- c(5,25,50,75,95)
    Qlabel <- c("5%","25%","50%","75%","95%")
    #  plot manifold along with marker and mesh points and crossing lines
    harmmat <- -fda::eval.fd(arclengthfine, harmvarmxfd)
    Qvec_al <- arclengthfine[ceiling(nfine*c(0.05, 0.25, 0.50, 0.75, 0.95))]
    Qharmmat <- -fda::eval.fd(Qvec_al, harmvarmxfd)
    if (nharm == 2) {
      df1 <- data.frame(harmmat)
      df2 <- data.frame(Qharmmat)
      pcaplt <- ggplot2::ggplot(df1, ggplot2::aes(harmmat[,1],harmmat[,2])) +
                ggplot2::geom_point(size=2) +
                ggplot2::geom_point(data=df2, ggplot2::aes(Qharmmat[,1],Qharmmat[,2], size=2)) +
                xlab('Rotated Component 1') +
                ylab('Rotated Component 2') +
                annotate("text", x=Qharmmat[1,1]*dodge, y=Qharmmat[1,2], label=Qlabel[1]) +
                annotate("text", x=Qharmmat[2,1]*dodge, y=Qharmmat[2,2], label=Qlabel[2]) +
                annotate("text", x=Qharmmat[3,1]*dodge, y=Qharmmat[3,2], label=Qlabel[3]) +
                annotate("text", x=Qharmmat[4,1]*dodge, y=Qharmmat[4,2], label=Qlabel[4]) +
                annotate("text", x=Qharmmat[5,1]*dodge, y=Qharmmat[5,2], label=Qlabel[5]) +
                theme(legend.position = "none",
                      axis.title=element_text(size=16,face="bold"))
      if (!is.null(titlestr))
      {
        pcaplt <- pcaplt + labs(title=titlestr)
      }
      
      print(pcaplt)
      
    } else {
      rgl::open3d()
      rgl::points3d( harmmat[,1],        harmmat[,2],  harmmat[,3], col = rainbow(1000), size=5) 
      rgl::points3d(Qharmmat[,1],       Qharmmat[,2], Qharmmat[,3], color="black", size=8)
      rgl::texts3d( Qharmmat[,1]*dodge, Qharmmat[,2], Qharmmat[,3], texts = Qlabel)
      rgl::axes3d()
      rgl::aspect3d(1,1,1)
      rgl::title3d(xlab="Rotated Component 1",ylab="Rotated Component 2",zlab="Rotated Component 3")
      pcaplt <- NULL
    }
    return(list(pcaplt=pcaplt, harmvarmxfd=harmvarmxfd, varpropvarmx=varpropvarmx))
  }
  
  
}