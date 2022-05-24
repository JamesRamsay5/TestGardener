Hfuns.plot <- function(theta, WfdList, U, plotindex=1) {
  evalarg  <- seq(0,100,len=51)
  Hval     <- Hfun(theta, WfdList, U)
  Result   <- DHfun(theta, WfdList, U)
  DHval    <- Result$DH
  D2Hval   <- Result$D2H
  linesize <- 1
  nindex   <- length(plotindex)
  plot_list <- list()
  for (j in 1:nindex) {
    indexj <- plotindex[j]
    Umatj  <- matrix(1,51,1) %*% U[indexj,]
    thetaj <- theta[indexj]
    Hj       <-  Hfun(evalarg, WfdList, Umatj)
    DHResult <- DHfun(evalarg, WfdList, Umatj)
    D2Hj <- DHResult$D2H
    #  plot function H
    df <- data.frame(x=evalarg,  y=Hj)
    p1 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  Hj)) +
      ggplot2::geom_line(size = linesize, color='blue') +
      ggplot2::geom_vline(xintercept = thetaj, size = linesize, color='blue', linetype = 2) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression(H(theta))) +
      ggplot2::labs(title=paste("Examinee",indexj,", theta =",round(thetaj, 2)))
    #  plot second derivative D2H
    df <- data.frame(x=evalarg,  y=D2Hj)
    p2 <- ggplot2::ggplot(df, ggplot2::aes(evalarg,  D2Hj)) +
      ggplot2::geom_line(size=linesize, color='blue') +
      ggplot2::geom_vline(xintercept = thetaj, color='blue', size=linesize, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, color='blue', size=linesize, linetype = 2) +
      ggplot2::xlab(expression(paste("Score index ", theta))) +
      ggplot2::ylab(expression(D2H(theta))) +
      ggplot2::labs(paste("Second derivative =",round(D2Hj,4)))
    p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)
    print(p)
    plot_list[[j]] <- p
    if (nindex > 1)  
      readline(prompt = paste("theta", indexj, ". Press [enter] to continue"))
  }
}

