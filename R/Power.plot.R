Power.plot <- function(WfdList, Qvec, dataList, plotindex=1:n, height=0.5,  
                       value=0, saveplot=FALSE, 
                       ttlsz=NULL, axisttl=NULL, axistxt=NULL) {
  n          <- length(WfdList)
  indfine    <- seq(0,100,len=101)
  optionList <- dataList$optList 
  ItemTitle  <- optionList$itemLab
  titlestr   <- dataList$titlestr
  nfine      <- length(indfine)
  linesize   <- 1
  powermat   <- matrix(0,nfine,n)
  powermatQ  <- matrix(0,n,5)
  plot_list  <- list()
  for (i in plotindex) {
    WListi <- WfdList[[i]]
    Wfdi   <- WListi$Wfd
    Mi     <- WListi$M
    for (k in 1:Mi) {
      powermat[,i] <- powermat[,i] + WListi$DWmatfine[,k]^2
    }
    powermat[,i] <- sqrt(powermat[,i])
    
    DWmatQ <- eval.surp(Qvec,Wfdi,1)
    for (k in 1:k) {
      powermatQ[i,] <- powermatQ[i,] + DWmatQ[,k]^2
    }
    powermatQ[i,] <- sqrt(powermatQ[i,])
    
    # ggplot version
    
    default_size=16
    default_size1=12
    default_size2=10
    My_Theme <- ggplot2::theme(
      plot.title = ggplot2::element_text(size = ifelse(is.null(ttlsz),  default_size,ttlsz)),
      axis.title = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
      axis.text  = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)))
    
    df <- data.frame(value=powermat[,i], indfine=indfine)
    pp <- ggplot2::ggplot(df, ggplot2::aes(indfine, powermat[,i])) +
      ggplot2::geom_line(size=linesize, na.rm=TRUE) +
      ggplot2::xlim(c(0,100)) + ggplot2::ylim(c(0,height)) +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
      ggplot2::annotate("point", x = Qvec[1], y = powermatQ[i,1], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[2], y = powermatQ[i,2], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[3], y = powermatQ[i,3], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[4], y = powermatQ[i,4], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[5], y = powermatQ[i,5], 
                        colour = "red", size = 1.5) +
      ggplot2::labs(title = paste("Question",i))+
      ggplot2::ylab(paste('Question power (',Mi,'-bits)',sep='')) +
      ggplot2::xlab("Score index") +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    
    # ----------- Juan Li 2021-02-17 ---------
    if (!is.null(titlestr))
    {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total power = ", 
                        round(pracma::trapz(indfine, powermat[,i]),2), sep="")
      } else 
      {
        ttllab <- paste('Question ', i,': ',ItemTitle[i], ", total power = ", 
                        round(pracma::trapz(indfine, powermat[,i]),2), sep="") 
        
      }
    } else
    {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total power = ", 
                        round(pracma::trapz(indfine, powermat[,i]),2), sep="")
      } else 
      {
        ttllab <- paste('Question ', i, ': ',ItemTitle[i], ", total power = ", 
                        round(pracma::trapz(indfine, powermat[,i]),2), sep="") 
        
      }
    }
    
    # -----------------------------------------
    
    
    pp <- pp + ggplot2::labs(title = ttllab)
    
    plot_list[[i]] <- pp
    print(pp)
    if (length(plotindex) > 1)
        readline(prompt = paste('Question', i, ". Press [enter] to continue"))
  }
  
  #  If required, a pdf file is set up containing all of the plots.  
  #  The file name is the string titlestr combined with "-power.pdf".  
  
  if (saveplot) {
    pdf(paste(titlestr,i,'-power.pdf',sep=""))
    for (i in plotindex) {
      print(plot_list[[i]])
    }
    dev.off()
  }
  
}

