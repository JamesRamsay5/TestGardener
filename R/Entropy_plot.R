Entropy_plot <- function(scrfine, SfdList, Qvec, dataList, plotindex=1:n, 
                         plotrange=c(min(scrfine),max(scrfine)),   
                         height=1.0, value=0,
                         ttlsz=NULL, axisttl=NULL, axistxt=NULL) {
  
  #  Last modified 9 November 23 by Jim Ramsay
  
  n          <- length(SfdList)
  optionList <- dataList$optList 
  ItemTitle  <- optionList$itemLab
  titlestr   <- dataList$titlestr
  nfine      <- length(scrfine)
  linesize   <- 1
  
  for (i in plotindex) {
    SListi <- SfdList[[i]]
    Sfdi   <- SListi$Sfd
    Mi     <- SListi$M
    Pmatfinei <- SListi$Pmatfine
    Smatfinei <- SListi$Smatfine
    
    #  compute entropy vector for this item over fine mesh
    
    entropyvec   <- rep(0,nfine)
    for (m in 1:Mi) {
      entropyvec<- entropyvec + Pmatfinei[,m]*Smatfinei[,m]
    }
    
    #  compute entropy values for marker percentages
    
    entropyvecQ  <- rep(0,5)
    SvecQ <- rep(0,5)
    for (m in 1:Mi) {
      for (k in 1:5) {
        SvecQ[k] <- pracma::interp1(as.numeric(scrfine), 
                                    as.numeric(Smatfinei[,m]), 
                                    as.numeric(Qvec[k]))
      }
      PvecQ <- Mi^(-SvecQ)
      entropyvecQ <- entropyvecQ + PvecQ*SvecQ
    }
    
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
    
    df <- data.frame(value=entropyvec, scrfine=scrfine)
    pp <- ggplot2::ggplot(df, ggplot2::aes(scrfine, entropyvec)) +
      ggplot2::geom_line(size=linesize, na.rm=TRUE) +
      ggplot2::xlim(plotrange) + ggplot2::ylim(c(0,height)) +
      ggplot2::geom_vline(xintercept = Qvec, color="black", 
                          linetype = "dashed") +
      ggplot2::annotate("point", x = Qvec[1], y = entropyvecQ[1], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[2], y = entropyvecQ[2], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[3], y = entropyvecQ[3], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[4], y = entropyvecQ[4], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[5], y = entropyvecQ[5], 
                        colour = "red", size = 1.5) +
      ggplot2::labs(title = paste("Question",i))+
      ggplot2::ylab(paste('Question entropy (',Mi,'-bits)',sep='')) +
      ggplot2::xlab("Score index") +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    
    if (!is.null(titlestr))
    {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total entropy = ", 
                        round(pracma::trapz(scrfine, entropyvec),2), sep="")
      } else {
        ttllab <- paste('Question ', i,': ',ItemTitle[i], ", total entropy = ", 
                        round(pracma::trapz(scrfine, entropyvec),2), sep="") 
      }
    } else {
      if (is.null(ItemTitle)) {
        ttllab <- paste('Question ', i, ", total entropy = ", 
                        round(pracma::trapz(scrfine, entropyvec),2), sep="")
      } else {
        ttllab <- paste('Question ', i, ': ',ItemTitle[i], ", total entropy = ", 
                        round(pracma::trapz(scrfine, entropyvec),2), sep="") 
      }
    }
    
    # -----------------------------------------
    
    pp <- pp + ggplot2::labs(title = ttllab)
    
    print(pp)
    
    if (length(plotindex) > 1) {
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

