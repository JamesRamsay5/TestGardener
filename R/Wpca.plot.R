library(plotly)

Wpca.plot <- function(arclength, WfdList, Wdim=NULL, nharm=2, rotate=TRUE, 
                       titlestr=NULL) {
  
  #  Last modified 13 July 2023 by Juan Li
  
  #  set up matrices of fine mesh values
  
  nfine   <- 101
  indfine <- seq(0,arclength,len=nfine)
  
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
  
  Pmat_full  <- matrix(0,nfine, Wdim)
  Wmat_full  <- matrix(0,nfine, Wdim)
  DWmat_full <- matrix(0,nfine, Wdim)
  
  n <- length(WfdList)
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
  
  #  set up basis for (smoothing over arclength)
  
  Wnbasis <- 7
  Wnorder <- 4
  Wbasis  <- fda::create.bspline.basis(c(0,arclength), Wnbasis, Wnorder) 
  WfdPar  <- fda::fdPar(fd(matrix(0,Wnbasis,1),Wbasis))
  
   #  fine mesh of points along manifold of length arclength
  
  nfine <- dim(Wmat_full)[1]
  arclengthfine <- seq(0,arclength,len=nfine)
  
  #  smooth all Wdim surprisal curves to get best fitting fd curves
  
  Sfd <- fda::smooth.basis(arclengthfine, Wmat_full, WfdPar)$fd
  
  #  functional PCA of the fd versions of surprisal curves
  #  the output is nharm principal component functions
  
  pcaList <- fda::pca.fd(Sfd, nharm, WfdPar, FALSE)
  
  #  set up the unrotated harmonic functional data object
  
  harmfd <- pcaList$harmonics
  
  #  compute variance proportions for unrotated solution
  
  varprop  <- pcaList$varprop
  eigvals  <- pcaList$values
  rooteigs <- sqrt(eigvals)[1:nharm]
  
  #  carry out a varimax rotation if desired
  
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
    #  plot manifold along with marker and mesh points and crossing lines
    harmmat  <- -fda::eval.fd(arclengthfine, harmvarmxfd)
    Qvec_al  <-  arclengthfine[ceiling(nfine*c(0.05, 0.25, 0.50, 0.75, 0.95))]
    Qvec_pts <- -fda::eval.fd(Qvec_al, harmvarmxfd)
    if (nharm == 2) {
      df_harmmat <- as.data.frame(harmmat)
      fig <- plot_ly(df_harmmat, x = ~PC1, y = ~PC2, type = 'scatter', 
                     mode = 'lines',
                     opacity = 1, 
                     line = list(width = 6, color = "black", reverscale = FALSE))%>% 
        layout(title = titlestr,
               showlegend = FALSE)
      
      df_Q <- as.data.frame(Qvec_pts)
      df_Q$label <- Qlabel
      # for (i in 1:5) {
      #   fig <- fig %>% 
      #     add_trace(x = df_Q$PC1[i], y = df_Q$PC2[i], 
      #               type = "scatter", text = df_Q$label[i], mode = "text")
      # }
      
    } else {
      df_harmmat <- as.data.frame(harmmat)
      fig <- plot_ly(df_harmmat, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', 
                     mode = 'lines',
                     opacity = 1, 
                     line = list(width = 6, color = "black", reverscale = FALSE))%>% 
        layout(title = titlestr,
               showlegend = FALSE)
      # Hi Juan,  The problem in rendering the info curve must be somewhere in this
      # code for plotting the marker quantiles because the error message is repeated
      # 5 times.  It occurs for both nharm = 2 and 3.
      # But when I comment out the following df_Q lines, nothing at all is rendered.
      # Could it be a problem with using the plotly layout?
      # df_Q <- as.data.frame(Qvec_pts)
      # df_Q$label <- Qlabel
      # message here:
      # A line object has been specified, but lines is not in the mode
      # Adding lines to the mode...
      # for (i in 1:5) {
      #   fig <- fig %>% 
      #     add_trace(x = df_Q$PC1[i], y = df_Q$PC2[i], z = df_Q$PC3[i], 
      #               type = "scatter3d", text = df_Q$label[i], mode = "text")
      # }
    }
  } else {
    stop("The current Wpca.plot only works with 2 or 3 dimension.")
  }
  
  return(list(pcaplt=fig, harmvarmxfd=harmvarmxfd, varpropvarmx=varpropvarmx,
              harmmat=harmmat, Qvec_pts=Qvec_pts, Qvec_al=Qvec_al))
  
}