Item.plot <- function(WfdList, Qvec, dataList, plotindex=1:n, key=NULL,
                      titlestr=NULL, saveplt=FALSE,
                      ttlsz=NULL,axisttl=NULL,axistxt=NULL,lgdlab=NULL,
                      width=c(-0.2,0.2), height=1)
{
  # Last modified 16 April 2021 by Juan Li
    
  n <- length(WfdList)
  if (is.null(plotindex)) plotindex <- 1:n
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  optionList <- dataList$optList # Juan Li 2021-02-18
  plot_list  <- list()
  
  powermat   <- matrix(0,nfine,n)
  powermatQ  <- matrix(0,n,5)
  
  entropymat   <- matrix(0,nfine,n)
  entropymatQ  <- matrix(0,n,5)
  
  linesize   <- 1
  
  ItemTitle  <- optionList$itemLab
  
  for (i in plotindex)
  {
    print(i)
    WListi    <- WfdList[[i]]
    Wfdi      <- WListi$Wfd
    Mi        <- WListi$M # 2020-12-14
    logMi     <- log(Mi)
    Wfitfinei <- WListi$Wmatfine
    Pfitfinei <- WListi$Pmatfine
    DWfitfinei <- WListi$DWmatfine
    itemStri  <- optionList$itemLab[i]  # Juan Li 2021-02-18
    optStri   <- optionList$optLab[[i]] # Juan Li 2021-02-18
    
    # power plot  
    for (k in 1:Mi) {
      powermat[,i] <- powermat[,i] + DWfitfinei[,k]^2
    }
    powermat[,i] <- sqrt(powermat[,i])
    
    DWmatQ <- eval.surp(Qvec,Wfdi,1)
    for (k in 1:k) {
      powermatQ[i,] <- powermatQ[i,] + DWmatQ[,k]^2
    }
    powermatQ[i,] <- sqrt(powermatQ[i,])
    
    # Entropy plot
    for (k in 1:Mi) {
      entropymat[,i] <- entropymat[,i] + Pfitfinei[,k]*Wfitfinei[,k]
    }
    
    WmatQ <- eval.surp(Qvec,Wfdi)
    PmatQ <- Mi^(-WmatQ)
    for (k in 1:Mi) {
      entropymatQ[i,] <- entropymatQ[i,] + PmatQ[,k]*WmatQ[,k]
    }
    
    if (is.null(dataList$key))
    {
      keyi <- NULL
    } else
    {
      keyi <- dataList$key[i]
    }
    
    if (!is.null(itemStri))
    {
      ttllab <- paste(titlestr,' ',i,': ',itemStri,sep="") # Juan Li 2021-02-18
    } else
    {
      ttllab <- paste(titlestr,': ','Question ', i,sep="")
    }
    
    # ----------- option label ---------
    if (!is.null(optStri)) 
    {
      optionVec <- optStri
    } else
    {
      optionVec <- NULL
    }
    # -----------------------------------------
    if (is.null(ItemTitle)) {
      ttllab <- paste('Question ', i, sep="")
    } else 
    {
      ttllab <- paste('Question ', i,': ',ItemTitle[i], sep="") 
    }
    
    # Probablity plot
    pp <-  plotP(Mi, binctr=NULL, indfine,
                Pbini=NULL, Pfitfinei, Qvec, keyi, FALSE, TRUE, 
                FALSE, FALSE, ttllab, optVec=optionVec) +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    
    # Sensitivity plot
    sp <- plotDW(Mi, indfine, DWfitfinei, Qvec, keyi, TRUE, FALSE, c(-0.2,0.2), ttllab,optVec=optionVec,
                ttlsz,axisttl,axistxt,lgdlab) +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    
    # Power plot
    df <- data.frame(value=powermat[,i], indfine=indfine)
    pp2 <- ggplot2::ggplot(df, ggplot2::aes(indfine, powermat[,i])) +
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
    
    ttllab_p <- paste("total power = ", 
                      round(pracma::trapz(indfine, powermat[,i]),2), sep="")
    pp2 <- pp2 + ggplot2::labs(title = ttllab_p)
    
    # Entropy plot
    df <- data.frame(value=entropymat[,i], indfine=indfine)
    ep <- ggplot2::ggplot(df, ggplot2::aes(indfine, entropymat[,i])) +
      ggplot2::geom_line(size=linesize, na.rm=TRUE) +
      ggplot2::xlim(c(0,100)) + ggplot2::ylim(c(0,height)) +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
      ggplot2::annotate("point", x = Qvec[1], y = entropymatQ[i,1], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[2], y = entropymatQ[i,2], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[3], y = entropymatQ[i,3], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[4], y = entropymatQ[i,4], 
                        colour = "red", size = 1.5) +
      ggplot2::annotate("point", x = Qvec[5], y = entropymatQ[i,5], 
                        colour = "red", size = 1.5) +
      ggplot2::labs(title = paste("Question",i))+
      ggplot2::ylab(paste('Question entropy (',Mi,'-bits)',sep='')) +
      ggplot2::xlab("Score index") +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))

    ttllab_e <- paste("total entropy = ", 
                      round(pracma::trapz(indfine, entropymat[,i]),2), sep="")
    ep <- ep + ggplot2::labs(title = ttllab_e)
    #-----------------------------------
    
    p <- ggpubr::ggarrange(pp, sp, pp2, ep, ncol = 2, nrow = 2,  legend = "right", common.legend = TRUE)
    plot_list[[i]] <- p
    if (length(plotindex) > 1)
    {
      print(p)
      readline(prompt = paste('Question', i, ". Press [enter] to continue"))
    }
    
  } # end i
  
  if (saveplt) {
    pdf(paste(titlestr,'-ICC.pdf',sep=""))
    for (i in plotindex) {
      print(plot_list[[i]])
    }
    dev.off()
  }
  
}
