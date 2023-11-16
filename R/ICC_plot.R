ICC_plot <- function(scrfine, SfdList, dataList, Qvec, 
                     binctr=NULL, data_point=FALSE, ci=FALSE, 
                     plotType="S", Srng=c(0,5), DSrng=c(-0.2, 0.2),
                     plotindex=1:n, 
                     titlestr=NULL, itemscopevec=rep(0,length(plotindex)),  
                     plotTitle=TRUE, autoplot=FALSE, plotMissing=TRUE, 
                     plotrange=c(min(scrfine),max(scrfine)), shaderange=NULL,  
                     ttlsz=NULL, axisttl=NULL, axistxt=NULL, 
                     lgdlab=NULL, lgdpos="bottom") {
  
  # plot a list of ICC plots
  
  # scrfine    Vector of length 101 containing plotting points
  # SfdList    List of Sbinsmth results
  # dataList   List of make.dataList
  # Qvec       vector of quantile positions
  
  # binctr     vector of bin centers, for plotting data points and/or 
  #            confidence interval
  # data_point default as FALSE; set to TRUE for plotting data points
  # ci         default as FALSE; set to TRUE for plotting confidence limits 
  # of probability and surprisal
  
  # plotType   Type(s) of plot, default as "P" for probability, 
  # can also be "S"  for surprisal, 
  #             "DS" for sensitivity, and any combination of the three
  
  # Srng       make ylim of all suprisal plots the same
  # DSrng      make ylim of all sensitivity plots the same
  
  # plotindex  a vector of item indices; set if uses want to focus on 
  # specific items 
  
  # titlestr   plot title
  # itemscopevec   vector of scope values
  # plotTitle  indicator of showing the plot title, default as TRUE
  
  # autoplot   plot all items in a batch
  
  # plotMissing Determine if plot the extra option for 
  # missing/spoiled responses.
  
  # plotrange  a vector of length 2; set if users want to focus on 
  # a specific score range 
  # shaderange a list of length 2 vector(s); set if users want to gray out 
  # specific score range(s)
  
  # ttlsz      font size of plot title
  # axisttl    font size of axis titles
  # axistxt    font size of axis labels
  # lgdlab     font size of legend labels
  # lgdpos     legend position, could be set as "None" to remove the legend
  
  #  --------------------------------------------------------------------------
  #               The arguments with defaults for debugging purposes
  #  ------------------------------------------------------------------------
  
  # binctr=NULL
  # data_point=FALSE
  # ci=FALSE
  # plotType="S"
  # Srng=c(0,5)
  # DSrng=c(-0.2, 0.2)
  # plotindex=1:n
  # titlestr=NULL
  # itemscopevec=rep(0,length(plotindex))
  # plotTitle=TRUE
  # autoplot=FALSE
  # plotMissing=TRUE
  # plotrange=c(min(scrfine),max(scrfine))
  # shaderange=NULL
  # ttlsz=NULL
  # axisttl=NULL
  # axistxt=NULL
  # lgdlab=NULL
  # lgdpos="bottom"
  
  # Last modified 1 November 2023 by Jim Ramsay
  
  n <- length(SfdList)

  #  --------------------------------------------------------------------------
  #                  Define background shading in plots
  #  --------------------------------------------------------------------------
  
  if (!is.null(shaderange)) {
    if (is.list(shaderange)) {
      shadeValid <- rep(TRUE, length(shaderange))
      for (ishade in 1:length(shaderange)) {
        rng <- shaderange[[ishade]]
        if (rng[1] >= rng[2]) {
          shadeValid[ishade] <- FALSE
        } else if (rng[1] < plotrange[1]) {
          shaderange[[ishade]][1] = plotrange[1]
        } else if (rng[2] > plotrange[2]) {
          shaderange[[ishade]][2] = plotrange[2]
        }
      }
      if (sum(shadeValid) > 0) {
        shaderangeV <- shaderange[shadeValid] 
      } else {
        shaderangeV <- NULL
      }
    } else {
      stop("shaderange is not a list.")
    }
  } else {
    shaderangeV <- NULL
  }

  #  --------------------------------------------------------------------------
  #               Loop through items to be plotted
  #  --------------------------------------------------------------------------
  
  nplot <- length(plotindex)
  plotlist <- list()
  
  #  ---------------------  autoplot && nplot > 1  --------------------
  
  if (autoplot && nplot > 1) {
    #  ------------------  plot all items in a batch. -----------------
    for (iplot in plotindex) {
      scopei <- itemscopevec[iplot]
      # print("autoplot")
      p <- plotCore(iplot, scrfine, SfdList, dataList, Qvec, 
                  binctr, data_point, ci, plotType, Srng, DSrng, 
                  titlestr, scopei, plotTitle, plotMissing, 
                  plotrange, shaderangeV, ttlsz, axisttl, axistxt,
                  lgdlab, lgdpos) 
      print(p)
      plotlist[[iplot]] <- p
    }
  } 
  
  #  ---------------------  length(plotindex) == 1  --------------------

    if (length(plotindex) == 1) {
      #  ------------------  plot all items in a batch. -----------------
        iplot <- plotindex
        # print(" length(plotindex) == 1")
        # print(plotindex)
        # print(data_point)
        # print(itemscopevec)
        scopei <- itemscopevec[iplot]
        p <- plotCore(iplot, scrfine, SfdList, dataList, Qvec, 
                      binctr, data_point, ci, plotType, Srng, DSrng, 
                      titlestr, scopei, plotTitle, plotMissing, 
                      plotrange, shaderangeV, ttlsz, axisttl, axistxt,
                      lgdlab, lgdpos) 
        print(p)
        plotlist[[1]] <- p
  } 
  
  #  ---------------------  !autoplot && nplot > 1  --------------------

    if (!autoplot && nplot > 1) {
    #  -------------  Usual case of user control of plotting  -----------
    # print("usual case")
    nplot <- length(plotindex)
    jplot <- 1
    iplot <- plotindex[jplot]
    while (jplot >= 1 & jplot <= nplot) {
      scopei <- itemscopevec[iplot]
      p <- plotCore(iplot, scrfine, SfdList, dataList, Qvec, 
                    binctr, data_point, ci, plotType, Srng, DSrng, 
                    titlestr, scopei, plotTitle, plotMissing, 
                    plotrange, shaderangeV, ttlsz, axisttl, axistxt,
                    lgdlab, lgdpos)  
      print(p)
      plotlist[[iplot]] <- p
      
      #  -------------------  request next plot. ----------------------
      
      if (length(plotindex) > 1) {
        n1 <- readline(prompt=
                         paste("Press [N/n/enter] to next item;", 
                               "[P/p] to previous item;", 
                               "or position in plotindex: "))
        if (n1 == "N" | n1 == "n" | n1 == "") {
          jplot <- jplot + 1
          iplot <- plotindex[jplot]
        } else {
          if (n1 == "P" | n1 == "p") {
            jplot <- jplot - 1
            iplot <- plotindex[jplot]
          } else {
            n1 <- suppressWarnings(as.integer(n1))
            if (!is.na(n1) & (n1 >= 1 & n1 <= nplot)) {
              jplot <- n1
              iplot <- plotindex[jplot]
            } else {
              break
            } 
          }
        }
      } else {
        break
      }
      gc() # to  free up memory from the previous run
    } 
    
  } # end plotting
  
  return(plotlist)
  
}

#. ----------------------------------------------------------------------------

plotCore <- function(iplot, scrfine, SfdList, dataList, Qvec, 
                     binctr, data_point, ci, plotType, Srng, DSrng, 
                     titlestr, scopei, plotTitle, plotMissing, 
                     plotrange, shaderange, ttlsz, axisttl, axistxt,
                     lgdlab, lgdpos) {
  
  # Last modified 13 October 2023 by Jim Ramsay
  
  #  obtain option labels
  
  # print("in plotCore")
  # print(data_point)
  optionList <- dataList$optList 
  
  #  obtain SfdList information for this plot
  
  SListi    <- SfdList[[iplot]]
  itemStri  <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
  optStri   <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
  Mi         <- SListi$M # 2020-12-14
  
  # obtain curve values at mesh points
  
  Sfitfinei  <- SListi$Smatfine
  DSfitfinei <- SListi$DSmatfine
  Pfitfinei  <- SListi$Pmatfine
  if (!plotMissing) {
    Mi         <- Mi - 1
    Sfitfinei  <- SListi$Smatfine[,1:Mi]
    DSfitfinei <- SListi$DSmatfine[,1:Mi]
    Pfitfinei  <- SListi$Pmatfine[,1:Mi]
  }
  
  # data for plotting data points
  
  if (data_point & !is.null(binctr)) {
    Pbini      <- SListi$Pbin
    Sbini      <- SListi$Sbin
    
    if (!plotMissing) {
      Pbini      <- SListi$Pbin[,1:Mi]
      Sbini      <- SListi$Sbin[,1:Mi]
    }
  } else {
    Pbini    <- NULL
    Sbini    <- NULL
  }
  
  # data for plotting confidence interval
  
  if (ci & !is.null(binctr)) {
    PStdErr    <- SListi$PStdErr   
    SStdErr    <- SListi$SStdErr
    
    if (!plotMissing)
    {
      PStdErr    <- SListi$PStdErr[,1:Mi]
      SStdErr    <- SListi$SStdErr[,1:Mi]
    }
  } else {
    PStdErr    <- NULL
    SStdErr    <- NULL
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
        if (scopei != 0) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr,         ' ',iplot,': ',itemStri,' ',"scope ",
                          round(scopei,1),sep="")
        } else {
          # print("no scope")
          # itemscopevec is missing
          ttllab <- paste(titlestr,' ',iplot,': ',itemStri,' ',"scope ",
                          round(scopei,1),sep="")
        }
      } else {
        # print("no itemstr")
        if (!(scopei == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr,         ' ',iplot,': ',"scope ",
                          round(scopei,1),sep="")
        } else {
          # print("no scope")
          # itemscopevec is missing
          ttllab <- paste(titlestr,' ',iplot, sep="")
        }
      }
    } else 
      if (!is.null(titlestr)) {
        # a scale title is supplied
        # print("titlestr")
        if (!(scopei == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste(titlestr," item ",iplot,": ","scope ",
                          round(scopei,1),sep="")
        } else {
          # print("no scope")
          # itemscopevec is missing
          ttllab <- paste(titlestr," item ",iplot,)
        }
      } else {
        # print("no titlestr")
        # print(scopei)
        if (!(scopei == 0)) {
          # print("scope")
          # a vector of scope value is supplied
          ttllab <- paste("item",iplot,': ',"scope ",
                          round(scopei,1))
        } else {
          # print("no scope")
          # itemscopevec is missing
          ttllab <- paste("item",iplot)
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
      # print("probability")
      ylabel <- "Proportion/Probability"
      pList[[itype]] <- plotICC(Mi, scrfine, Pfitfinei, Qvec, keyi,
                                binctr, Pbini, PStdErr, c(0,1), 0.5, 
                                ttllab, xlabel, ylabel, optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "S") {
      # surprisal
      # print("invoke plotICC with surprisal")
      ylabel <- paste("Surprisal (",Mi,"-bits)",sep="")
      pList[[itype]] <- plotICC(Mi, scrfine, Sfitfinei, Qvec, keyi,
                                binctr, Sbini, SStdErr, Srng, 0, 
                                ttllab, xlabel, ylabel, optionVec,
                                plotrange, shaderange, 
                                ttlsz, axisttl, axistxt,
                                lgdlab, lgdpos)
    } else if (plotType[itype] == "DS") {
      # sensitivity
      ylabel <- "Sensitivity"
      pList[[itype]] <- plotICC(Mi, scrfine, DSfitfinei, Qvec, keyi,
                                binctr, NULL, NULL, DSrng, 0, 
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

#. ----------------------------------------------------------------------------

plotICC   <- function(Mi, scrfine, fitfinei, Qvec, keyi, 
                      binctr, bin1, StdErr, 
                      range, intercept, 
                      ttllab, xlab, ylab, optVec,
                      plotrange, shaderange, 
                      ttlsz, axisttl, axistxt,
                      lgdlab, lgdpos)
{
  # plot a single ICC plot
  # print("in plotICC")
  # Last modified 27 April 2023 by Juan Li
  value    <- 0
  variable <- 0
  if (!is.null(binctr)) nbin <- length(binctr)
  linesize <- 1
  ptssize  <- 2
  ind  <- 1:Mi 
  if (is.null(keyi)) indS <- NULL  else indS  <- ind[ind!=keyi] 
  
  # ----------- Juan Li 2021-02-17 ---------
  if (is.null(optVec)) 
  {
    optVec <- rep("",Mi)
    for(m in 1:Mi)
      optVec[m] <- m #paste("Option: ",m, sep="")
  }
  
  ymin <- NULL
  ymax <- NULL
  ymin_r <- NULL
  ymax_r <- NULL
  
  if (!is.null(keyi)) {
    # print("key is not null")
    # print(keyi)
    #  ---------------------------------------------------------------------
    #       Plot probability curve(s) for multiple choice test items
    #  ---------------------------------------------------------------------
    
    # plot the curve of right option
    dffine_r <- data.frame(scrfine=scrfine,
                           fitfinei_r= fitfinei[,keyi])
    pp <- ggplot2::ggplot(
      data=dffine_r, ggplot2::aes(scrfine,fitfinei[,keyi])) +
      ggplot2::geom_line(ggplot2::aes(lty=optVec[keyi]),
                         color="blue", linewidth=linesize*2, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = intercept,  color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")
    
    # plot the curve of wrong options
    dffine_w <- data.frame(scrfine=scrfine, fitfinei_w= fitfinei[,indS])
    names(dffine_w)[2:ncol(dffine_w)] <- optVec[indS]
    dffine_w <- tidyr::gather(dffine_w, key = "variable", value = "value", -scrfine)
    
    pp <- pp +
      ggplot2::geom_line(data = dffine_w, ggplot2::aes(scrfine,value,color=variable), 
                         linewidth=linesize, na.rm = TRUE) +
      ggplot2::scale_colour_hue(name="Srong")
    
  } else {
    # print("key null")
    #  ---------------------------------------------------------------------
    #             Plot probability curves for scale items
    #  ---------------------------------------------------------------------
    
    dffine <- data.frame(scrfine=scrfine, fitfinei = fitfinei[,ind]) 
    names(dffine)[2:ncol(dffine)]<- optVec[ind]
    dffine <-   tidyr::gather(dffine, key = "variable", value = "value", -scrfine)
    
    pp <- ggplot2::ggplot(dffine, ggplot2::aes(scrfine,value,color=variable)) +
      ggplot2::geom_line(linewidth=linesize, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = intercept,  color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
      ggplot2::scale_colour_hue()
  }
  
  # plot data points (optional)
  if (!is.null(binctr) & !is.null(bin1)) # plot observed proportions
  {
    if (!is.null(keyi)) {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for multiple choice test items
      #  ---------------------------------------------------------------------
      
      # plot data points of right option
      dfpts_r <- data.frame(binctr=binctr,
                            bin1_r= bin1[,keyi])
      
      pp <- pp + ggplot2::geom_point(
        data=dfpts_r, ggplot2::aes(binctr,bin1[,keyi]), shape = 21, fill = "blue", size = ptssize, 
        na.rm = TRUE)
      
      # plot data points of wrong options
      dfpts_w <- data.frame(binctr=binctr,
                            bin1_w= bin1[,indS]) 
      names(dfpts_w)[2:ncol(dfpts_w)]<- optVec[indS]
      dfpts_w <-   tidyr::gather(dfpts_w, key = "variable", value = "value", -binctr)
      
      pp <- pp +
        ggplot2::geom_point(data=dfpts_w, ggplot2::aes(binctr,value,fill=variable), 
                            shape = 21, size = ptssize, na.rm = TRUE)+
        ggplot2::scale_fill_hue(guide = "none")
    } else {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for scale items
      #  ---------------------------------------------------------------------
      
      dfpts <- data.frame(binctr=binctr,
                          bin1= bin1[,ind]) 
      names(dfpts)[2:ncol(dfpts)]<- optVec[ind]
      dfpts <-   tidyr::gather(dfpts, key = "variable", value = "value", -binctr)
      
      pp <- pp +
        ggplot2::geom_point(data=dfpts, ggplot2::aes(binctr,value,fill=variable),
                            shape = 21, size = ptssize, na.rm = TRUE)+
        ggplot2::scale_fill_hue(guide = "none")
    }
  }
  
  # plot confidence interval (optional)
  if (!is.null(StdErr) & !is.null(binctr))  
  {
    binfit <- matrix(0,nbin,Mi)
    for (m in 1:Mi) {
      binfit[,m] <- pracma::interp1(as.numeric(scrfine), as.numeric(fitfinei[,m]), 
                                    as.numeric(binctr))
    }
    
    if (!is.null(keyi)) {
      
      #  ---------------------------------------------------------------------
      #       Plot confidence intervals for multiple choice test items
      #  ---------------------------------------------------------------------
      
      # ci of thr right option
      ci_r <- data.frame(binctr = binctr,
                         ymin_r = binfit[,keyi]-2*StdErr[,keyi],
                         ymax_r = binfit[,keyi]+2*StdErr[,keyi])
      
      pp <- pp + ggplot2::geom_ribbon(data = ci_r, ggplot2::aes(ymin = ymin_r, ymax = ymax_r, x = binctr), 
                                      inherit.aes = FALSE, fill = "blue", alpha = 0.5)
      
      # ci of wrong options
      cimin_w <- data.frame(binctr = binctr,
                            ymin = binfit[,indS]-2*StdErr[,indS]) 
      names(cimin_w)[2:ncol(cimin_w)]<- optVec[indS]
      cimin_w <-   tidyr::gather(cimin_w, key = "variable", value = "ymin", -binctr)
      
      cimax_w <- data.frame(binctr = binctr,
                            ymax   = binfit[,indS]+2*StdErr[,indS]) 
      names(cimax_w)[2:ncol(cimax_w)]<- optVec[indS]
      cimax_w <-   tidyr::gather(cimax_w, key = "variable", value = "ymax", -binctr)
      
      ci_w    <- cimin_w %>% dplyr::left_join(cimax_w, by = c("binctr", "variable"))
      
      pp <- pp +
        ggplot2::geom_ribbon(data=ci_w, ggplot2::aes(ymin = ymin, ymax = ymax, x = binctr, fill=variable), 
                             inherit.aes = FALSE, alpha = 0.5)+
        ggplot2::scale_fill_hue(guide = "none")
      
    } else {
      
      #  ---------------------------------------------------------------------
      #       Plot confidence intervals for scale items
      #  ---------------------------------------------------------------------
      
      cimin <- data.frame(binctr = binctr,
                          ymin   = binfit[,ind]-2*StdErr[,ind]) 
      names(cimin)[2:ncol(cimin)]<- optVec[ind]
      cimin <-   tidyr::gather(cimin, key = "variable", value = "ymin", -binctr)
      
      cimax <- data.frame(binctr = binctr,
                          ymax   = binfit[,ind]+2*StdErr[,ind]) 
      names(cimax)[2:ncol(cimax)]<- optVec[ind]
      cimax <-   tidyr::gather(cimax, key = "variable", value = "ymax", -binctr)
      
      ci    <- cimin %>% dplyr::left_join(cimax, by = c("binctr", "variable"))
      
      pp <- pp +
        ggplot2::geom_ribbon(data=ci, ggplot2::aes(ymin = ymin, ymax = ymax, x = binctr, fill=variable), 
                             inherit.aes = FALSE, alpha = 0.5)+
        ggplot2::scale_fill_hue(guide = "none")
    }
  }
  
  if (!is.null(shaderange))
  {
    rectAlpha <- 0.9
    for (ishade in 1:length(shaderange))
    {
      pp <- pp + ggplot2::annotate("rect", xmin=shaderange[[ishade]][1], xmax=shaderange[[ishade]][2], 
                                   ymin=range[1], ymax=range[2],
                                   alpha = rectAlpha)
    }
  }
  
  pp <- pp +
    ggplot2::ylim(range[1], range[2]) +
    ggplot2::xlim(plotrange[1], plotrange[2]) +
    ggplot2::labs(x = xlab, 
                  y = ylab)
  
  if (!is.null(ttllab)) {
    pp <- pp +
      ggplot2::labs(title = ttllab)
  }
  
  default_size =16
  default_size1=12
  default_size2=10
  
  pp <- pp + ggplot2::theme(legend.position = lgdpos,
                            legend.title=ggplot2::element_blank(),
                            plot.title  = ggplot2::element_text(size = ifelse(is.null(ttlsz),  default_size,ttlsz)),
                            axis.title  = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
                            axis.text   = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)),
                            legend.text = ggplot2::element_text(size = ifelse(is.null(lgdlab), default_size2,lgdlab))) # remove legend title
  
  return(pp)
}


