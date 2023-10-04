plotCore   <- function(iplot, scrfine, WfdList, dataList, Qvec, 
                       binctr, data_point, ci, plotType, Wrng, DWrng, 
                       titlestr, scopeveci, plotTitle, plotMissing, 
                       plotrange, shaderange, ttlsz, axisttl, axistxt,
                       lgdlab, lgdpos) {
  
  # Last modified 2 October 2023 by Jim Ramsay
  
  #  obtain option labels
  
  optionList <- dataList$optList 
  
  #  obtain WfdList information for this plot
  
  WListi    <- WfdList[[iplot]]
  itemStri  <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
  optStri   <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
  Mi         <- WListi$M # 2020-12-14
  
  # obtain curve values at mesh points
  
  Wfitfinei  <- WListi$Wmatfine
  DWfitfinei <- WListi$DWmatfine
  Pfitfinei  <- WListi$Pmatfine
  if (!plotMissing) {
    Mi         <- Mi - 1
    Wfitfinei  <- WListi$Wmatfine[,1:Mi]
    DWfitfinei <- WListi$DWmatfine[,1:Mi]
    Pfitfinei  <- WListi$Pmatfine[,1:Mi]
  }
  
  # data for plotting data points
  
  if (data_point & !is.null(binctr)) {
    Pbini      <- WListi$Pbin
    Wbini      <- WListi$Wbin
    
    if (!plotMissing) {
      Pbini      <- WListi$Pbin[,1:Mi]
      Wbini      <- WListi$Wbin[,1:Mi]
    }
  } else {
    Pbini    <- NULL
    Wbini    <- NULL
  }
  
  # data for plotting confidence interval
  
  if (ci & !is.null(binctr)) {
    PStdErr    <- WListi$PStdErr   
    WStdErr    <- WListi$WStdErr
    
    if (!plotMissing)
    {
      PStdErr    <- WListi$PStdErr[,1:Mi]
      WStdErr    <- WListi$WStdErr[,1:Mi]
    }
  } else {
    PStdErr    <- NULL
    WStdErr    <- NULL
  }
  
  # vector for answer key
  
  if (is.null(dataList$key)) {
    keyi <- NULL
  } else {
    keyi <- dataList$key[iplot]
  }
  
  # ------------------  assemble the title for this item  ----------------------
  
  if (plotTitle) {
    # print("plot title")
    # assemble the title ttllab for this item
    if (!is.null(itemStri)) {
      # an item title string is supplied
      # print("itemStr")
      if (!is.null(titlestr)) {
        # a scale title is supplied
        # print("titlestr")
        if (!(scopeveci == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr,         ' ',iplot,': ',itemStri,' ',"scope ",
                          round(scopeveci,1),sep="")
        } else {
          # print("no scope")
          # scopevec is missing
          ttllab <- paste(dataList$titlestr,' ',iplot,': ',itemStri,' ',"scope ",
                          round(scopeveci,1),sep="")
        }
      } else {
        # print("no itemstr")
        if (!(scopeveci == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr,         ' ',iplot,': ',"scope ",
                          round(scopeveci,1),sep="")
        } else {
          # print("no scope")
          # scopevec is missing
          ttllab <- paste(titlestr,' ',iplot, sep="")
        }
      }
    } else {
      if (!is.null(titlestr)) {
        # a scale title is supplied
        # print("titlestr")
        if (!(scopeveci == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr," item ",iplot,": ","scope ",
                          round(scopeveci,1),sep="")
        } else {
          # print("no scope")
          # scopevec is missing
          ttllab <- paste(titlestr," item ",iplot,)
        }
      } else {
        # print("no titlestr")
        # print(scopeveci)
        if (!(scopeveci == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste("item",iplot,': ',"scope ",
                          round(scopeveci,1),sep="")
        } else {
          # print("no scope")
          # scopevec is missing
          ttllab <- paste("item",iplot)
        }
      }
    }
  } else {
    # don't plot a title
    ttllab <- NULL
  }
  
  # --------------  option string. -------------
  
  if (!is.null(optStri)) {
    optionVec <- optStri
  } else {
    optionVec <- NULL
  }
  
  # ---------------  construct abscissa label xlabel. -------------
  
  if (max(scrfine) == 100) {
    xlabel <- "Score Index" 
  } else {
    xlabel <- paste("Information Metric (",Mi," bits)", sep="")    
  }
  
  #  ---------------- Call plotICC to make the plot. ---------------------------
  
  nplotType  <- length(plotType)
  pList <- list()
  for (itype in 1:nplotType) {
    if (plotType[itype] == "P") {
      # probability
      ylabel <- "Proportion/Probability"
      pList[[itype]] <- plotICC(Mi, scrfine, Pfitfinei, Qvec, keyi,
                                binctr, Pbini, PStdErr, c(0,1), 0.5, 
                                ttllab, xlabel, ylabel, optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "S" || plotType[itype] == "W") {
      # surprisal
      ylabel <- paste("Surprisal (",Mi,"-bits)",sep="")
      pList[[itype]] <- plotICC(Mi, scrfine, Wfitfinei, Qvec, keyi,
                                binctr, Wbini, WStdErr, Wrng, 0, 
                                ttllab, xlabel, ylabel, optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "DS" || plotType[itype] == "DW") {
      # sensitivity
      ylabel <- "Sensitivity"
      pList[[itype]] <- plotICC(Mi, scrfine, DWfitfinei, Qvec, keyi,
                                binctr, NULL, NULL, DWrng, 0, 
                                ttllab, xlabel, ylabel, optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else {
      stop("Can't recognize the plot type.")
    }
  }
  
  if (nplotType > 1) {
    p <- ggpubr::ggarrange(plotlist = pList, ncol = 1)
  } else {
    p <- pList[[1]]
  }
  
  return(p)
}