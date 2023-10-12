Wpca.plot <- function(harmvarmxfd, nharm = 2, titlestr = NULL) {
  
  #  Last modified 10 October 2023 by Jim Ramsay
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  #  display the approximate test manifold curve if required
  #  The coordinates of a point are the nharm harmonic values at that point
  
  #  plot the first two or three harmonics using ggplot2 or plotly
  
  if (nharm == 2 || nharm == 3) {
    Qvec       <- c(5,25,50,75,95)
    harmmat    <- -fda::eval.fd(indfine, harmvarmxfd)
    Qvec_pts   <- -fda::eval.fd(Qvec,    harmvarmxfd)
    df_harmmat <- as.data.frame(harmmat)
    df_Q       <- as.data.frame(Qvec_pts)
    df_Q$label <- c("5%","25%","50%","75%","95%")
    lineList <- list(width = 6, color = "black", reverscale = FALSE)
    
    #. --------------------------  set up the plot ------------------------
    if (nharm == 2) {
      #   -----------------------  two dimensions -----------------------
      plot_ly(df_harmmat, x = ~PC1, y = ~PC2, type = 'scatter', 
        mode = 'lines', opacity = 1, line = lineList) %>% 
        layout(title = titlestr, showlegend = FALSE)  %>% 
        add_trace(x = df_Q$PC1[1], 
                  y = df_Q$PC2[1], 
                  type = "scatter", 
                  text = df_Q$label[1], 
                  mode = "text+lines")  %>%
      add_trace(x = df_Q$PC1[2], 
                y = df_Q$PC2[2], 
                type = "scatter", 
                text = df_Q$label[2], 
                mode = "text+lines")  %>%
      add_trace(x = df_Q$PC1[3], 
                y = df_Q$PC2[3], 
                type = "scatter", 
                text = df_Q$label[3], 
                mode = "text+lines")  %>%
      add_trace(x = df_Q$PC1[4], 
                y = df_Q$PC2[4], 
                type = "scatter", 
                text = df_Q$label[4], 
                mode = "text+lines")  %>%
      add_trace(x = df_Q$PC1[5], 
                y = df_Q$PC2[5], 
                type = "scatter", 
                text = df_Q$label[5], 
                mode = "text+lines")
    } else {
      #   -----------------------  three dimensions -----------------------
      plot_ly(df_harmmat, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', 
        mode = 'lines', opacity = 1, line = lineList) %>% 
        layout(title = titlestr, showlegend = FALSE)  %>%
        add_trace(x = df_Q$PC1[1], 
                  y = df_Q$PC2[1], 
                  z = df_Q$PC3[1], 
                  type = "scatter3d", 
                  text = df_Q$label[1], 
                  mode = "text+lines")  %>%
        add_trace(x = df_Q$PC1[2], 
                  y = df_Q$PC2[2], 
                  z = df_Q$PC3[2], 
                  type = "scatter3d", 
                  text = df_Q$label[2], 
                  mode = "text+lines")   %>%
        add_trace(x = df_Q$PC1[3], 
                  y = df_Q$PC2[3], 
                  z = df_Q$PC3[3], 
                  type = "scatter3d", 
                  text = df_Q$label[3], 
                  mode = "text+lines")  %>%
        add_trace(x = df_Q$PC1[4], 
                  y = df_Q$PC2[4], 
                  z = df_Q$PC3[4], 
                  type = "scatter3d", 
                  text = df_Q$label[4], 
                  mode = "text+lines")  %>%
        add_trace(x = df_Q$PC1[5], 
                  y = df_Q$PC2[5], 
                  z = df_Q$PC3[5], 
                  type = "scatter3d", 
                  text = df_Q$label[5], 
                  mode = "text+lines")
    }
  } else {
    stop("The current Wpca.plot only works with 2 or 3 dimension.")
  }
  
}