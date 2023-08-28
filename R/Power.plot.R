Power.plot <- function(scrfine, WfdList, Qvec, dataList, plotindex=1:n, 
                       plotrange=c(min(scrfine),max(scrfine)),   
                       height=0.5, value=0,
                       ttlsz=NULL, axisttl=NULL, axistxt=NULL) {
  
  #  Last modified 20 May 2022 by Jim Ramsay
  
  n          <- length(WfdList)
  optionList <- dataList$optList 
  ItemTitle  <- optionList$itemLab
  titlestr   <- dataList$titlestr
  nfine      <- length(scrfine)
  linesize   <- 1
  powermatQ  <- matrix(0,    5,n)
  
  for (i in plotindex) {
    WListi     <- WfdList[[i]]
    Wfdi       <- WListi$Wfd
    Mi         <- WListi$M
    DWmatfinei <- WListi$DWmatfine
    
    #  compute power curve values over fine mesh
    
    powervec <- rep(0,nfine)
    for (m in 1:Mi) {
      powervec <- powervec + DWmatfinei[,m]^2
    }
    powervec <- sqrt(powervec)
    
    powervecQ <- rep(0,5)
    DWvecQ <- rep(0,5)
    for (m in 1:Mi) {
      for (k in 1:5) {
        DWvecQ[k] <- pracma::interp1(as.numeric(scrfine), 
                                     as.numeric(DWmatfinei[,m]), 
                                     as.numeric(Qvec[k]))
      }
      powervecQ <- powervecQ + DWvecQ^2
    }
    powervecQ <- sqrt(powervecQ)
    
    # ggplot version
    
    default_size =16
    default_size1=12
    default_size2=10
    My_Theme <- ggplot2::theme(
      plot.title = 
        ggplot2::element_text(size = 
                                ifelse(is.null(ttlsz),  default_size,ttlsz)),
      axis.title = 
        ggplot2::element_text(size = 
                                ifelse(is.null(axisttl),default_size1,axisttl)),
      axis.text  = 
        ggplot2::element_text(size = 
                                ifelse(is.null(axistxt),default_size2,axistxt)))
    
    df <- data.frame(value=powervec, scrfine=scrfine)
    pp <- ggplot2::ggplot(df, ggplot2::aes(scrfine, powervec)) +
      ggplot2::geom_line(size=linesize, na.rm=TRUE) +
      ggplot2::xlim(plotrange) + ggplot2::ylim(c(0,height)) +
      ggplot2::geom_vline(xintercept = Qvec, color="black", 
                          linetype = "dashed") +
      ggplot2::annotate("point", x = Qvec[1], y = powervecQ[1], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[2], y = powervecQ[2], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[3], y = powervecQ[3], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[4], y = powervecQ[4], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[5], y = powervecQ[5], 
                        colour = "red", size = 1.5) +
      ggplot2::labs(title = paste("Question",i))+
      ggplot2::ylab(paste('Question power (',Mi,'-bits)',sep='')) +
      ggplot2::xlab("Score index") +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    
    if (!is.null(titlestr))
    {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total power = ", 
                        round(pracma::trapz(scrfine, powervec),2), sep="")
      } else {
        ttllab <- paste('Question ', i,': ',ItemTitle[i], ", total power = ", 
                        round(pracma::trapz(scrfine, powervec),2), sep="") 
      }
    } else {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total power = ", 
                        round(pracma::trapz(scrfine, powervec),2), sep="")
      } else {
        ttllab <- paste('Question ', i, ': ',ItemTitle[i], ", total power = ", 
                        round(pracma::trapz(scrfine, powervec),2), sep="") 
      }
    }
    
    # -----------------------------------------
    
    pp <- pp + ggplot2::labs(title = ttllab)
    
    print(pp)
    
    if (length(plotindex) > 1)
    {
      n1<-readline(prompt=
            paste("Press [N/n/enter] to next item;",
                  "[P/p] to previous item; or item index to that item: "))
      if (n1 == "N" | n1 == "n" | n1 == "") {
        i = i + 1
      } else if (n1 == "P" | n1 == "p") {
        i = i - 1
      } else {
        n1 <- suppressWarnings(as.integer(n1))
        if (!is.na(n1) & (n1 >= plotindex[1] & n1 <= 
                          plotindex[length(plotindex)])) {
          i <- n1
        } else {
          break
        }
      }
    }
    gc() # to  free up memory from the previous run
  }
}