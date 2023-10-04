Wpca.plot <- function(WfdList, Wdim=NULL, 
                      nharm=2, rotate=TRUE, titlestr=NULL) {
  
  #  Last modified 27 September 2023 by Jim Ramsay
  
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
  
  #  display the approximate test manifold curve if required
  #  The coordinates of a point are the nharm harmonic values at that point
  
  #  plot the first two or three harmonics using ggplot2 or plotly
  
  if (nharm == 2 || nharm == 3) {
    pind   <- c(5,25,50,75,95)
    Qlabel <- c("5%","25%","50%","75%","95%")
    harmmat    <- -fda::eval.fd(indfine, harmvarmxfd)
    Qvec <- pind
    Qvec_pts   <- -fda::eval.fd(Qvec, harmvarmxfd)
    df_harmmat <- as.data.frame(harmmat)
    df_Q       <- as.data.frame(Qvec_pts)
    df_Q$label <- Qlabel
    if (nharm == 2) {
      fig <- plot_ly(df_harmmat, x = ~PC1, y = ~PC2, type = 'scatter', 
                     mode = 'lines',
                     opacity = 1, 
                     line = list(width = 6, color = "black", reverscale = FALSE))%>% 
        layout(title = titlestr,
               showlegend = FALSE)
      
      for (i in 1:5) {
        fig <- fig %>% 
          add_trace(x = df_Q$PC1[i], y = df_Q$PC2[i], 
                    type = "scatter", text = df_Q$label[i], mode = "text+lines")
      }
      
    } else {
      fig <- plot_ly(df_harmmat, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', 
                     mode = 'lines',
                     opacity = 1, 
                     line = list(width = 6, color = "black", reverscale = FALSE) ) %>% 
        layout(title = titlestr,
               showlegend = FALSE)
      
      for (i in 1:5) {
        fig <- fig %>% 
          add_trace(x = df_Q$PC1[i], y = df_Q$PC2[i], z = df_Q$PC3[i], text = df_Q$label[i], 
                    type = "scatter3d", mode = "text+lines")
      }
    }
  } else {
    stop("The current Wpca.plot only works with 2 or 3 dimension.")
  }
  
  return(list(pcaplt=fig, harmvarmxfd=harmvarmxfd, varpropvarmx=varpropvarmx,
              harmmat=harmmat, Qvec_pts=Qvec_pts, Qvec=Qvec))
  
}