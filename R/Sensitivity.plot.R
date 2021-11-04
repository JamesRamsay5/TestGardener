Sensitivity.plot <- function(WfdList, Qvec, dataList, plotindex=1:n, key=NULL,
                             titlestr=NULL, saveplot=FALSE, width=c(-0.2,0.2),
                             ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
{
	
  #  Last modified 28 March 2021 by Jim Ramsay
  
  n <- length(WfdList)
  if (is.null(plotindex))
  {
    plotindex <- 1:n
  }
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine) 
  optionList <- dataList$optList # Juan Li 2021-02-18
  plot_list = list()
  
  for (iplot in plotindex)
  {
    WListi     <- WfdList[[iplot]]
    Wfdi       <- WListi$Wfd
    Mi         <- WListi$M # 2020-12-14
    DWfitfinei <- eval.surp(indfine, Wfdi, 1)
    itemStri   <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
    optStri    <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
    
    if (is.null(key))
    {
      keyi <- NULL
    } else
    {
      keyi <- key[iplot]
    }
    
    if (!is.null(itemStri))
    {
      ttllab <- paste(titlestr,' ',iplot,': ',itemStri,sep="") # Juan Li 2021-02-18
    } else
    {
      ttllab <- paste(titlestr,': ','Question ', iplot,sep="")
    }
    
    # ----------- Juan Li 2021-02-17 ---------
    if (!is.null(optStri)) 
    {
      optionVec <- optStri
    } else
    {
      optionVec <- NULL
    }
    # -----------------------------------------
    
    alltype   <- TRUE # logical, if true, plot right and wrong P curves
    mssplot   <- FALSE # logical, if true, plot spoiled responses
    
    p <- plotDW(Mi, indfine, DWfitfinei, Qvec, keyi, alltype, mssplot, width, 
                ttllab, optVec=optionVec,
                ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
    
    plot_list[[iplot]] <- p
    print(p)
    if (length(plotindex) > 1)
        readline(prompt = paste('Question', iplot, ". Press [enter] to continue"))
  } # end iplot
  
  if (saveplot) {
    pdf(paste(titlestr,'-sensitivity.pdf',sep=""))
    for (i in plotindex) {
      print(plot_list[[i]])
    }
    dev.off()
  }
}

#  ----------------------------------------------------------------------------

plotDW   <- function(Mi, indfine, DWfitfinei, Qvec, keyi=NULL,  
                     alltype=TRUE, mssplot=FALSE, DWrng, ttllab, optVec=NULL,
                     ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
{
  value    <- 0
  variable <- 0
  linesize <- 1
  ptssize  <- 2
  ind  <- 1:Mi 
  if (is.null(keyi)) indW <- NULL  else indW  <- ind[ind!=keyi] 
  
  # ----------- Juan Li 2021-02-17 ---------
  if (is.null(optVec)) 
  {
    optVec <- rep("",Mi)
    for(m in 1:Mi)
      optVec[m] <- paste("Option: ",m, sep="")
  }
  # -----------------------------------------
  default_size =16
  default_size1=12
  default_size2=10
  My_Theme = ggplot2::theme(
    plot.title  = ggplot2::element_text(size = ifelse(is.null(ttlsz),  default_size,ttlsz)),
    axis.title  = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
    axis.text   = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)),
    legend.text = ggplot2::element_text(size = ifelse(is.null(lgdlab), default_size2,lgdlab)))
  
  # Plot ICC curve(s)
  if (!is.null(keyi))
  {
    
    #  ---------------------------------------------------------------------
    #       Plot sensitivity curve(s) for multiple choice test items
    #  ---------------------------------------------------------------------
    
    # plot the curve of right option
    dffine_r <- data.frame(indfine=indfine,
                           DWfitfinei_r=DWfitfinei[,keyi])
    
    dwp <- ggplot2::ggplot(data=dffine_r, ggplot2::aes(indfine,DWfitfinei[,keyi])) +
           ggplot2::geom_line(ggplot2::aes(lty=optVec[keyi]),
                              color="blue", size=linesize*2, na.rm = TRUE) +
           ggplot2::scale_linetype('Right') +
           ggplot2::geom_hline(yintercept = 0,    color="black", linetype = "dashed") +
           ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")
    
    if (alltype)# also plot curves of wrong options
    {
      dffine_w <- data.frame(indfine=indfine,
                             DWfitfinei_w= DWfitfinei[,indW])
      names(dffine_w)[2:ncol(dffine_w)]<- optVec[indW]
      dffine_w <- tidyr::gather(dffine_w, key = "variable", value = "value", -indfine)
      
      dwp <- dwp +
        ggplot2::geom_line(data = dffine_w, ggplot2::aes(indfine,value,color=variable),
                           size=linesize, na.rm = TRUE) +
        ggplot2::scale_colour_hue(name="Wrong")
    }
  } else {

    #  ---------------------------------------------------------------------
    #             Plot sensitivity curves for scale items
    #  ---------------------------------------------------------------------

    dffine <- data.frame(indfine=indfine, DWfitfinei_r = DWfitfinei[,ind])
    names(dffine)[2:ncol(dffine)]<- optVec[ind]
    dffine <-  tidyr::gather(dffine, key = "variable", value = "value", -indfine)
    
    dwp <- ggplot2::ggplot(dffine, ggplot2::aes(indfine,value,color=variable)) +
           ggplot2::geom_line(size=linesize, na.rm = TRUE) +
           ggplot2::geom_hline(yintercept = 0,    color="black", linetype = "dashed") +
           ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
           ggplot2::scale_colour_hue(name="Wrong")
      
  }
  dwp <- dwp + ggplot2::labs(title = ttllab)
  
  if (!is.null(DWrng))
  {
    dwp <- dwp + ggplot2::ylim(DWrng)
  }
  
  dwp <- dwp +
    ggplot2::xlim(0,100) +
    ggplot2::xlab("Score index") +
    ggplot2::ylab("Sensitivity") +
    ggplot2::theme(axis.title  =ggplot2::element_text(size=16,face="bold")) +
    ggplot2::theme(legend.title=ggplot2::element_blank()) # remove legend title
  
  return(dwp)
}

