Wbinsmth.plot <- function(scrfine, plotindex=1:n, 
                          plotrange=c(min(scrfine),max(scrfine)), 
                          binctr, Qvec, dataList, WfdList, twoplot=TRUE, 
                          ptsplot=TRUE, alltype=TRUE, landscape=FALSE,  
                          Wrng=c(0,5), cnfplot=FALSE, ttlsz=NULL, 
                          axisttl=NULL, axistxt=NULL, lgdlab=NULL) {
  
  # Last modified 18 May 2022 by Jim Ramsay
  
  n <- length(WfdList)
  if (is.null(plotindex)) plotindex <- 1:n
  
  optionList <- dataList$optList # Juan Li 2021-02-18
  plot_list  <- list()
  
  for (iplot in plotindex)
  {
    WListi    <- WfdList[[iplot]]
    Wfdi      <- WListi$Wfd
    Mi        <- WListi$M # 2020-12-14
    logMi     <- log(Mi)
    Pbini     <- WListi$Pbin
    Wbini     <- WListi$Wbin
    Pfitfinei <- WListi$Pmatfine
    Wfitfinei <- WListi$Wmatfine
    itemStri  <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
    optStri   <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
    
    if (is.null(dataList$key))
    {
      keyi <- NULL
    } else
    {
      keyi <- dataList$key[iplot]
    }
    
    if (!is.null(itemStri))
    {
      ttllab <- paste(dataList$titlestr,' ',iplot,': ',itemStri,sep="") 
    } else
    {
      ttllab <- paste(dataList$titlestr,': ','Question ', iplot,sep="")
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
    
    if (!twoplot) {  
     # One-panel plot
      p <-  plotP(Mi, binctr, scrfine,
                  Pbini, Pfitfinei, Qvec, keyi, ptsplot, alltype, 
                  twoplot, landscape, optVec=optionVec)
    } else {  
     # Two-panel plot
      pp <- plotP(Mi, binctr, scrfine,
                  Pbini, Pfitfinei, Qvec, keyi, ptsplot, alltype, 
                  twoplot, landscape, optVec=optionVec)
      wp <- plotW(Mi, binctr, scrfine,
                  Wbini, Wfitfinei, Qvec, keyi, ptsplot, alltype, 
                  twoplot, landscape, optVec=optionVec, Wrng)

      if (landscape) # plots are side by side
      {
        p <- ggpubr::ggarrange(pp, wp, ncol = 2, nrow = 1, legend = "right", 
                               common.legend = TRUE)
      } else
      {
        p <- ggpubr::ggarrange(pp, wp, ncol = 1, nrow = 2,  legend = "right", 
                               common.legend = TRUE)
      }
    }
    
    print(p)
    
    plot_list[[iplot]] <- p
    if (length(plotindex) > 1)
      readline(prompt = paste('Question', iplot, ". Press [enter] to continue"))
      
  } # end iplot
  
  # if (saveplot) {
  #   pdf(paste(dataList$titlestr,'-ICC.pdf',sep=""))
  #   for (i in plotindex) {
  #     print(plot_list[[i]])
  #   }
  #   dev.off()
  # }
  
}

#  ----------------------------------------------------------------------------

plotP   <- function(Mi, binctr, scrfine, Pbini, Pfitfinei, 
                    Qvec, keyi=NULL, ptsplot=TRUE, alltype=TRUE, 
                    twoplot=FALSE, landscape=FALSE, optVec=NULL)
{
  value    <- 0
  variable <- 0
  nbin     <- length(binctr)
  linesize <- 1
  ptssize  <- 2
  ind      <- 1:Mi 
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
    plot.title  = ggplot2::element_text(size = default_size),
    axis.title  = ggplot2::element_text(size = default_size1),
    axis.text   = ggplot2::element_text(size = default_size2),
    legend.text = ggplot2::element_text(size = default_size2))
  
  if (!is.null(keyi)) {
    
    #  ---------------------------------------------------------------------
    #       Plot probability curve(s) for multiple choice test items
    #  ---------------------------------------------------------------------
    
    # plot the curve of right option
    dffine_r <- data.frame(scrfine=scrfine,
                           Pfitfinei_r= Pfitfinei[,keyi])
    pp <- ggplot2::ggplot(
      data=dffine_r, ggplot2::aes(scrfine,Pfitfinei[,keyi])) +
      ggplot2::geom_line(ggplot2::aes(lty=optVec[keyi]),
                color="blue", size=linesize*2, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = 0.5,  color="black", 
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", 
                                              linetype = "dashed")
    
    if (alltype)
    {
      #  alltype is TRUE, plot all curves
      dffine_w <- data.frame(scrfine=scrfine, Pfitfinei_w= Pfitfinei[,indW])
      names(dffine_w)[2:ncol(dffine_w)] <- optVec[indW]
      dffine_w <- tidyr::gather(dffine_w, key = "variable", 
                                value = "value", -scrfine)
      
      pp <- pp +
        ggplot2::geom_line(data = dffine_w, 
                           ggplot2::aes(scrfine,value,color=variable), 
                                        size=linesize, na.rm = TRUE) +
                           ggplot2::scale_colour_hue(name="Wrong")
                           # palette=ggplot2::scales::hue_pal(direction = -1)
    }
    
  } else {
    
    #  ---------------------------------------------------------------------
    #             Plot probability curves for scale items
    #  ---------------------------------------------------------------------
    
    dffine <- data.frame(scrfine=scrfine, Pfitfinei = Pfitfinei[,ind]) 
    names(dffine)[2:ncol(dffine)]<- optVec[ind]
    dffine <- tidyr::gather(dffine, key = "variable", value = "value", -scrfine)
    
    pp <- 
      ggplot2::ggplot(dffine, ggplot2::aes(scrfine,value,color=variable)) +
      ggplot2::geom_line(size=linesize, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = 0.5,  color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
      ggplot2::scale_colour_hue()
      # palette=ggplot::scales::hue_pal(direction = -1)
    
  }
  
  # Plot data points (optional)
  
  if (ptsplot) # plot observed proportions
  {
    if (!is.null(keyi)) {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for multiple choice test items
      #  ---------------------------------------------------------------------
      
      # plot data points of right option
      dfpts_r <- data.frame(binctr=binctr,
                            Pbini_r= Pbini[,keyi])
      
      pp <- pp + ggplot2::geom_point(data=dfpts_r, 
                 ggplot2::aes(binctr,Pbini[,keyi]), shape = 21, 
                          fill = "blue", size = ptssize, 
        na.rm = TRUE)
      
      if (alltype) # plot data points of wrong options
      {
        dfpts_w <- data.frame(binctr=binctr,
                              Pbini_w= Pbini[,indW]) 
        names(dfpts_w)[2:ncol(dfpts_w)]<- optVec[indW]
        dfpts_w <-   tidyr::gather(dfpts_w, key = "variable",  
                                   value = "value", -binctr)
        
        pp <- pp +
          ggplot2::geom_point(data=dfpts_w, 
                              ggplot2::aes(binctr,value,fill=variable), 
                                           shape = 21, size = ptssize, 
                                           na.rm = TRUE) +
          ggplot2::scale_fill_hue(guide = "none")
          # palette=sggplot2::cales::hue_pal(direction = -1)
      }
    } else {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for scale items
      #  ---------------------------------------------------------------------
      
      dfpts <- data.frame(binctr=binctr,
                          Pbini= Pbini[,ind]) 
      names(dfpts)[2:ncol(dfpts)]<- optVec[ind]
      dfpts <-   tidyr::gather(dfpts, key = "variable", 
                               value = "value", -binctr)
      
      pp <- pp +
        ggplot2::geom_point(data=dfpts, 
                            ggplot2::aes(binctr,value,fill=variable),
                                         shape = 21, size = ptssize, 
                                         na.rm = TRUE) +
        ggplot2::scale_fill_hue(guide = "none")
        # palette=ggplot2::scales::hue_pal(direction = -1)
    }
  }
  
  # pp <- pp + ggplot2::labs(title = ttllab)
  pp <- pp +
    ggplot2::ylim(0,1) +
    ggplot2::xlim(0,max(scrfine)) +
    ggplot2::xlab("Score Index") +
    ggplot2::ylab("Proportion/Probability")
  
  if (twoplot) {
    if (landscape) {
      pp <- pp +
        ggplot2::xlab("Score index") +
        ggplot2::theme(axis.title=ggplot2::element_text(size=16, face="bold"),
                       aspect.ratio=1)
    } else {
      pp <- pp +
        ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
    }
  } else {
    pp <- pp +
      ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
  }
  
  # remove legend title
  pp <- pp + ggplot2::theme(legend.title=ggplot2::element_blank()) 
  pp <- pp + My_Theme
  return(pp)
}

#  -----------------------------------------------------------------

plotW   <- function(Mi, binctr, scrfine,
                    Wbini, Wfitfinei, Qvec,
                    keyi=NULL, ptsplot=TRUE, alltype=TRUE, 
                    twoplot=FALSE, landscape=FALSE, optVec=NULL, Wrng)
{
  value    <- 0
  variable <- 0
  nbin <- length(binctr)
  linesize <- 1
  ptssize  <- 2
  ind  <- 1:Mi 
  if (is.null(keyi)) indW <- NULL
  else indW  <- ind[ind!=keyi] # indices of wrong options
  
  # ----------- Juan Li 2021-02-17 ---------
  if (is.null(optVec)) 
  {
    optVec <- rep("",Mi)
    for(m in 1:Mi)
      optVec[m] <- paste("Option: ",m, sep="")
  }
  # -----------------------------------------
  default_size=16
  default_size1=12
  default_size2=10
  My_Theme = ggplot2::theme(
    plot.title  = ggplot2::element_text(size = default_size ),
    axis.title  = ggplot2::element_text(size = default_size1),
    axis.text   = ggplot2::element_text(size = default_size2),
    legend.text = ggplot2::element_text(size = default_size2))
  
  # Plot ICC curve(s)
  if (!is.null(keyi))
  {
    # plot the curve of right option
    dffine_r <- data.frame(scrfine=scrfine,
                           Wfitfinei_r= Wfitfinei[,keyi])
    
    wp <- ggplot2::ggplot(data=dffine_r, 
                          ggplot2::aes(scrfine,Wfitfinei[,keyi])) +
      ggplot2::geom_line(ggplot2::aes(lty=optVec[keyi]),
                         color="blue", size=linesize*2, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = 0,    color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")
    
    if (alltype)# also plot curves of wrong options
    {
      dffine_w <- data.frame(scrfine=scrfine,
                             Wfitfinei_w= Wfitfinei[,indW])
      names(dffine_w)[2:ncol(dffine_w)] <- optVec[indW]
      dffine_w <- tidyr::gather(dffine_w, key = "variable", 
                                value = "value", -scrfine)
      
      wp <- wp +
        ggplot2::geom_line(data = dffine_w, 
                           ggplot2::aes(scrfine,value,color=variable), 
                                        size=linesize, na.rm = TRUE) +
        ggplot2::scale_colour_hue(name="Wrong")
        # palette=scales::hue_pal(direction = -1)
    }
  } else
  {
    dffine <- data.frame(scrfine=scrfine, Wfitfinei=Wfitfinei[,ind])
    names(dffine)[2:ncol(dffine)] <- optVec[ind]
    dffine <- tidyr::gather(dffine, key ="variable", 
                            value ="value", -scrfine)
    
    wp <- ggplot2::ggplot(dffine, 
                          ggplot2::aes(scrfine,value,color=variable)) +
          ggplot2::geom_line(size=linesize, na.rm = TRUE) +
          ggplot2::geom_hline(yintercept = 0, color="black", 
                              linetype = "dashed") +
          ggplot2::geom_vline(xintercept = Qvec, color="black", 
                              linetype = "dashed") +
          ggplot2::scale_colour_hue()
          # palette=ggplot2::scales::hue_pal(direction = -1)
  }
  
  # Plot data points (optional)
  if (ptsplot) # plot observed proportions
  {
    if (!is.null(keyi))
    {
      # plot data points of right option
      dfpts_r <- data.frame(binctr=binctr,
                          Wbini_r= Wbini[,keyi])
      
      wp <- wp + 
            ggplot2::geom_point(data=dfpts_r, 
                                ggplot2::aes(binctr,Wbini[,keyi]), 
                                     shape = 21, fill = "blue", 
                                     size = ptssize, na.rm = TRUE)
      
      if (alltype)# plot right and wrong P curves
      {
        dfpts_w <- data.frame(binctr=binctr,
                              Wbini_w= Wbini[,indW]) 
        names(dfpts_w)[2:ncol(dfpts_w)]<- optVec[indW]
        dfpts_w <-   tidyr::gather(dfpts_w, key = "variable", 
                                   value = "value", -binctr)
        
        wp <- wp +
          ggplot2::geom_point(data=dfpts_w, 
                              ggplot2::aes(binctr,value,fill=variable), 
                                           shape = 21, size = ptssize, 
                                           na.rm = TRUE) +
          ggplot2::scale_fill_hue(guide = "none")
          # palette=scales::hue_pal(direction = -1)
      }
    } else
    {
      dfpts <- data.frame(binctr=binctr,
                          Wbini= Wbini[,ind])
      names(dfpts)[2:ncol(dfpts)]<- optVec[ind]
      dfpts <-  tidyr::gather(dfpts, key ="variable", 
                              value ="value", -binctr)
      
      wp <- wp +
        ggplot2::geom_point(data=dfpts, 
                            ggplot2::aes(binctr,value,fill=variable), 
                                         shape = 21, size = ptssize, 
                                         na.rm = TRUE) +
        ggplot2::scale_fill_hue(guide = "none")
        # palette=scales::hue_pal(direction = -1)
    }
  }
  
  if (!is.null(Wrng))
  {
    wp <- wp + ggplot2::ylim(Wrng)
  }
  
  ystr <- paste("Surprisal (",Mi,"-bits)",sep="")
  
  wp <- wp +
    ggplot2::xlim(min(scrfine),max(scrfine)) +
    ggplot2::xlab("Score index") +
    ggplot2::ylab(ystr)+
    ggplot2::theme(axis.title=ggplot2::element_text(size=12,face="bold"))
  
  if (twoplot)
  {
    if (landscape)
    {
      wp <- wp +
        ggplot2::theme(aspect.ratio=1)
    }
  }
  # remove legend title
  wp <- wp + ggplot2::theme(legend.title=ggplot2::element_blank()) 
  wp <- wp + My_Theme
  return(wp)
}

