scoreDensity <- function(scrvec, scrrng=c(0,100), ndensbasis=15, 
                         ttlstr=NULL, pltmax=0) {

# Last modified 9 November 2023 by Jim Ramsay

  #  define a fine mesh of score values for plotting
  
  nfine   <- 101
  scrfine <- seq(scrrng[1],scrrng[2],len=101)
  
  #  set up the bspline basis object
  
  logdensbasis <- fda::create.bspline.basis(scrrng, ndensbasis)

  # call index_distn to get the smooth density function
  
  tdList       <- index_distn(scrvec, logdensbasis, nfine)
  
  denscdf      <- tdList$cdffine
  logdensfd    <- tdList$logdensfd
  logdensvec   <- fda::eval.fd(scrfine, logdensfd)
  densfine     <- exp(logdensvec)/tdList$C
  pvec         <- c(0.05,0.25,0.50,0.75,0.95)
  denscdfi     <- unique(denscdf)
  scrfinei     <- seq(scrrng[1],scrrng[2],length.out=length(denscdfi))
  Qvec         <- pracma::interp1(as.numeric(denscdfi), 
                                  as.numeric(scrfinei), as.numeric(pvec))
  
  #  plot using indexPlot below
  
  dens.plot <- index_density.plot(scrvec, scrrng, scrfine, densfine, Qvec, 
                                  ttlstr, TRUE, pltmax)
  
  print(dens.plot)
  
  return(dens.plot)
}

#  ----------------------------------------------------------------------------

index_density.plot <- function(scrvec, scrrng, scrfine, densfine, Qvec, 
                              ttlstr=NULL, labelwrd = TRUE, pltmax=0) {
  
  # Last modified 9 February 2021 by Jim Ramsay
  
  binstr   <-   floor(scrrng[1])
  binend   <- ceiling(scrrng[2])
  binwidth <- 1

  if (pltmax==0) pltmax <- max(densfine)
  if (labelwrd) {
    label_y  <- pltmax/20
  } else {
    label_y  <- 0
  }
  # print(label_y)
  
  # Plot proportions
  
  ..count.. <- NULL
  df1 <- data.frame(scrvec=scrvec)
  df2 <- data.frame(scrfine=scrfine, densfine=densfine)
  dens.plot <- ggplot2::ggplot(df1, ggplot2::aes(scrvec)) +
    ggplot2::geom_histogram(binwidth = binwidth, color="black", fill=NA, 
    ggplot2::aes(y=(..count..)/sum(..count..)), na.rm = TRUE)
  dens.plot <- dens.plot +  
    ggplot2::geom_line(data=df2, ggplot2::aes(x=scrfine,y=densfine, color="red"), 
                       linewidth=2,na.rm = TRUE) +  
    ggplot2::geom_vline(xintercept = Qvec, linetype = 2) +
    ggplot2::ylim(0,pltmax) +
    ggplot2::xlab("score") +
    ggplot2::ylab("proportion") +
    ggplot2::theme(legend.position = "none",
                 axis.title=ggplot2::element_text(face="bold")) 
  dens.plot <- dens.plot +
    ggplot2::annotate("text", x=Qvec[1], y=label_y, label=" 5%") +
    ggplot2::annotate("text", x=Qvec[2], y=label_y, label="25%") +
    ggplot2::annotate("text", x=Qvec[3], y=label_y, label="50%") +
    ggplot2::annotate("text", x=Qvec[4], y=label_y, label="75%") +
    ggplot2::annotate("text", x=Qvec[5], y=label_y, label="95%") 
  if (!is.null(ttlstr)) dens.plot <- dens.plot + ggplot2::labs(title = ttlstr)
  
  return(dens.plot)
}

