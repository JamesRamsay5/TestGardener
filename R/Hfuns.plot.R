Hfuns.plot <- function(evalarg, theta, WfdList, U, plotindex=1) {
  #  Plot a person's fit function
  #  Arguments:
  #  evalarg   ... A vector containing values over which the plotting takes place.
  #  theta     ... A vector of score index values
  #  WfdList   ... A list vector containing descriptions of each item
  #  U         ... The N by n matrix of choice indices
  #  plotindex ... The the indices in theta of cases to be plotted
  
  #  Last modified 7 August 2023 by Jim Ramsay
  
  if (!is.matrix(evalarg)) evalarg <- matrix(evalarg, 101, 1)
  nevalarg <- length(evalarg)
  if (!is.matrix(U)) U <- t(as.matrix(U))
  linesize <- 1
  nplot   <- length(plotindex)
  plot_list <- list()
  for (j in 1:nplot) {
    indexj    <- plotindex[j]
    thetaj    <- theta[indexj]
    Umatj     <- as.numeric(U[indexj,])
    Hj        <- TestGardener::Hcurve(WfdList, Umatj)
    D2Hj      <- TestGardener::Hcurve(WfdList, Umatj, nderiv=2)
    D2Hj[101] <- D2Hj[100]
    #  plot function H
    df <- data.frame(x=evalarg,  y=Hj)
    p1 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  Hj)) +
          ggplot2::geom_line(linewidth = linesize, color='blue') +
          ggplot2::geom_vline(xintercept = thetaj, 
                              linewidth = linesize, color='blue', linetype = 2) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression(H(theta))) +
      ggplot2::labs(title=paste("Examinee",indexj,", theta =",round(thetaj, 2)))
    #  plot second derivative D2H
    df <- data.frame(x=evalarg,  y=D2Hj)
    p2 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  D2Hj)) +
      ggplot2::geom_line(linewidth=linesize, color='blue') +
      ggplot2::geom_vline(xintercept = thetaj, color='blue', 
                          linewidth=linesize, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, color='blue', 
                          linewidth=linesize, linetype = 2) +
      ggplot2::xlab(expression(paste("Score index ", theta))) +
      ggplot2::ylab(expression(D2H(theta))) +
      ggplot2::labs(paste("Second derivative =",round(D2Hj,4)))
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)
    print(p)
    plot_list[[j]] <- p
    if (nplot > 1)
      readline(prompt = paste("theta", indexj, ". Press [enter] to continue"))
  }
}

