plotCore   <- function(iplot, scrfine, WfdList, dataList, Qvec, 
                       binctr, data_point, ci,
                       plotType, 
                       Wrng, DWrng, 
                       titlestr, plotTitle, xlab, ylab,
                       plotMissing, 
                       plotrange, shaderange,  
                       ttlsz, axisttl, axistxt,
                       lgdlab, lgdpos) {
  
  # Last modified 10 July 2023 by Juan Li
  
  optionList <- dataList$optList 
  
  WListi    <- WfdList[[iplot]]
  itemStri  <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
  optStri   <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
  Mi         <- WListi$M # 2020-12-14
  
  # data for plotting the curves
  Wfitfinei  <- WListi$Wmatfine
  DWfitfinei <- WListi$DWmatfine
  Pfitfinei  <- WListi$Pmatfine
  if (!plotMissing)
  {
    Mi         <- Mi - 1
    Wfitfinei  <- WListi$Wmatfine[,1:Mi]
    DWfitfinei <- WListi$DWmatfine[,1:Mi]
    Pfitfinei  <- WListi$Pmatfine[,1:Mi]
  }
  
  # data for plotting data points
  if (data_point & !is.null(binctr)) 
  {
    Pbini      <- WListi$Pbin
    Wbini      <- WListi$Wbin
    
    if (!plotMissing)
    {
      Pbini      <- WListi$Pbin[,1:Mi]
      Wbini      <- WListi$Wbin[,1:Mi]
    }
  } else
  {
    Pbini    <- NULL
    Wbini    <- NULL
  }
  
  # data for plotting confidence interval
  if (ci & !is.null(binctr)) 
  {
    PStdErr    <- WListi$PStdErr   
    WStdErr    <- WListi$WStdErr
    
    if (!plotMissing)
    {
      PStdErr    <- WListi$PStdErr[,1:Mi]
      WStdErr    <- WListi$WStdErr[,1:Mi]
    }
  } else
  {
    PStdErr    <- NULL
    WStdErr    <- NULL
  }
  
  # vector of answer key
  if (is.null(dataList$key))
  {
    keyi <- NULL
  } else
  {
    keyi <- dataList$key[iplot]
  }
  
  # plot title for each item
  if (plotTitle) {
    if (!is.null(itemStri))
    {
      if (!is.null(titlestr))
      {
        ttllab <- paste(titlestr,' ',iplot,': ',itemStri,sep="")
      } else
      {
        ttllab <- paste(dataList$titlestr,' ',iplot,': ',itemStri,sep="")
      }
    } else
    {
      if (!is.null(titlestr))
      {
        ttllab <- paste(titlestr,' ', iplot,sep="")
      } else
      {
        ttllab <- paste(dataList$titlestr,' ', iplot,sep="")
      }
    }
  } else {
    ttllab <- NULL
  }
  
  # option string
  if (!is.null(optStri)) 
  {
    optionVec <- optStri
  } else
  {
    optionVec <- NULL
  }
  
  # -----------------------------------------
  nplot  <- length(plotType)
  pList <- list()
  
  # ylabel
  ylabel <- rep("", nplot)
  for (itype in 1:nplot) {
    if (plotType[itype] == "P") {
      ylabel[itype] = "Proportion/Probability"
    } else if (plotType[itype] == "W") {
      ylabel[itype] = paste("Surprisal (",Mi,"-bits)",sep="")
    } else {
      ylabel[itype] = "Sensitivity"
    }
  }
  if (length(ylab) == nplot) { # 2023-05-28,  & nplot > 1
    ylabel[1] = ylab
  }
  
  #  Call plotICC
  for (itype in 1:nplot) {
    if (plotType[itype] == "P") {
      # probability
      pList[[itype]] <- plotICC(Mi, scrfine, Pfitfinei, Qvec, keyi,
                                binctr, Pbini, PStdErr, 
                                c(0,1), 0.5, 
                                ttllab, xlab, ylabel[itype], optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "W") {
      # surprisal
      pList[[itype]] <- plotICC(Mi, scrfine, Wfitfinei, Qvec, keyi,
                                binctr, Wbini, WStdErr, 
                                Wrng, 0, 
                                ttllab, xlab, ylabel[itype], optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "DW") {
      # sensitivity
      pList[[itype]] <- plotICC(Mi, scrfine, DWfitfinei, Qvec, keyi,
                                binctr, NULL, NULL, 
                                DWrng, 0, 
                                ttllab, xlab, ylabel[itype], optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else {
      stop("Can't recognize the plot type.")
    }
  }
  
  if (nplot > 1) {
    p <- ggpubr::ggarrange(plotlist = pList, ncol = 1)
  } else {
    p <- pList[[1]]
  }
  
  return(p)
}