surprisal.curve <- function(ST, bdry) {
  
  # Display a three dimensional space curve within a 3d surprisal mesh
  
  #. Arguments:
  
  #  ST   ... An n by 3 matrix containing a sequence of 3 surprisal
  #           values within a 3d surprisal surface
  #  bdry ... outer limit of the surprisal mesh
  
  #  Last modified 28 February 2023
  
  #  Select log base M
  M <- 3
  #   define orthonormal map from dimension 3 to dimension 2
  Zmat <- zerobasis(M) 
  #  set up 101 points for plotting coordinate lines
  coordvalues <- seq(-bdry,bdry, length.out = 101)
  #  set up range a number of vector values
  n <- 21
  yvalues <- seq(-bdry, bdry, length.out = n)
  #   generate the planar mesh from the n y-value vector
  res <- pracma::meshgrid(yvalues, yvalues)
  B1 <- res$X
  B2 <- res$Y
  #  reshape into matrix with n^2 rows and two columns
  B <- cbind(as.vector(B1),as.vector(B2))
  #  convert B to surprisal and probability 3D values
  res <- surprisal.chart(B, Zmat, M)
  S <- res$S
  #  change surprisal values back into n by n matrices
  S1 <- matrix(S[,1],nrow = n)
  S2 <- matrix(S[,2],nrow = n)
  S3 <- matrix(S[,3],nrow = n)
  #  convert to 3D values
  res <- surprisal.chart(cbind(coordvalues, rep(0, 101)), Zmat, M)
  SX <- res$S
  res <- surprisal.chart(cbind(rep(0, 101), coordvalues), Zmat, M)
  SB <- res$S
  # trim plot
  trim = TRUE
  if (trim)
  {
    max <- max(SX)
    S1[S1>max] <- NA
    S2[S2>max] <- NA
    S3[S3>max] <- NA
  }
  
  figS <- plotly::plot_ly()
  for (i in seq(1,n)) {
    figS <- plotly::add_trace(figS, x = S1[,i], y = S2[,i], z = S3[,i],
                              visible = TRUE, type = 'scatter3d', 
                              mode = 'lines', line = list(color = "blue"))
    figS <- plotly::add_trace(figS, x = S1[,i], y = S3[,i], z = S2[,i],
                              visible = TRUE, type = 'scatter3d', 
                              mode = 'lines', line = list(color = "blue"))
  }
  figS <- plotly::add_trace(figS, x = SX[,1], y = SX[,2], z = SX[,3],
                            visible = TRUE, type = 'scatter3d', 
                            mode = 'lines', 
                            line = list(color = "black", width = 4))
  figS <- plotly::add_trace(figS, x = SB[,1], y = SB[,2], z = SB[,3],
                            visible = TRUE, type = 'scatter3d', 
                            mode = 'lines', 
                            line = list(color = "black", width = 4))
  axx <- list(nticks = 3, title = "Surprisal 1")
  axy <- list(nticks = 3, title = "Surprisal 2")
  axz <- list(nticks = 3, title = "Surprisal 3")
  figS <- figS %>% plotly::layout(showlegend = FALSE,
                   scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  #  plot 3D surprisal surface with trajectory
  #  ST <- eval.surp(seq(0,100, length.out = 1001), STfd$Wfd)
  # add curve to figS
  STE <- nrow(ST)
  figS <- plotly::add_trace(figS, x = ST[   ,1], y = ST[   ,2], z = ST[   ,3],
                            visible = TRUE, type = 'scatter3d', 
                            mode = 'lines', 
                            line   = list(color = "red", width = 4))
  #  plot curve origin
  figS <- plotly::add_trace(figS, x = ST[  1,1], y = ST[  1,2], z = ST[  1,3],
                            visible = TRUE, type = 'scatter3d', 
                            mode = 'markers',
                            marker = list(color = "red", size = 8, 
                                          symbol = 'circle'))
  #  plot curve termination
  figS <- plotly::add_trace(figS, x = ST[STE,1], y = ST[STE,2], z = ST[STE,3],
                            visible = TRUE, type = 'scatter3d', 
                            mode = 'markers',
                            marker = list(color = "red", size = 4, 
                                          symbol = 'x'))
  
  figS

  
}