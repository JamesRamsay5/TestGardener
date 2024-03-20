Ffuns_plot <- function(evalarg, index, SfdList, chcemat, plotindex=1) {
  #  Plot a person's fit function
  #  Arguments:
  #  evalarg   ... A vector containing values over which the plotting takes place.
  #  index     ... A vector of score index values
  #  SfdList   ... A list vector containing descriptions of each item
  #  chcemat   ... The N by n matrix of choice indices
  #  plotindex ... The the indices in index of cases to be plotted
  
  #  Last modified 9 November 2023 by Jim Ramsay
  
  if (!is.matrix(evalarg)) evalarg <- matrix(evalarg, 101, 1)
  nevalarg <- length(evalarg)
  if (!is.matrix(chcemat)) chcemat <- t(as.matrix(chcemat))
  linesize <- 1
  nplot     <- length(plotindex)
  plot_list <- list()
  for (j in 1:nplot) {
    plotj      <- plotindex[j]
    indexj     <- plotj
    chcematj   <- as.numeric(chcemat[indexj,])
    Fj         <- Fcurve(SfdList, chcematj)
    D2Fj       <- Fcurve(SfdList, chcematj, nderiv=2)
    D2Fj[101]  <- D2Fj[100]
    #  plot function F
    df <- data.frame(x=evalarg,  y=Fj)
    p1 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  Fj)) +
      ggplot2::geom_line(linewidth = linesize, color='blue') +
      ggplot2::geom_vline(xintercept = index[plotj], 
                          linewidth = linesize, color='blue', linetype = 2) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression(F(index))) +
      ggplot2::labs(title=paste("Examinee ",plotj,", index = ",round(index[plotj], 1),sep=""))
    #  plot second derivative D2F
    df <- data.frame(x=evalarg,  y=D2Fj)
    p2 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  D2Fj)) +
      ggplot2::geom_line(linewidth=linesize, color='blue') +
      ggplot2::geom_vline(xintercept = index[plotj], color='blue', 
                          linewidth=linesize, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, color='blue', 
                          linewidth=linesize, linetype = 2) +
      ggplot2::xlab(expression(Score(index))) +
      ggplot2::ylab(expression(D2F(index))) +
      ggplot2::labs(paste("Second derivative =",round(D2Fj,4)))
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)
    print(p)
    plot_list[[j]] <- p
    if (nplot > 1)
      readline(prompt = paste("Examinee", round(plotj,1), "  Press [enter] to continue"))
  }
}

