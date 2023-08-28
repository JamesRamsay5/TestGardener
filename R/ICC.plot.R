ICC.plot <- function(scrfine, WfdList, dataList, Qvec, 
                     binctr = NULL, data_point = FALSE, ci = FALSE, 
                     plotType = "P", 
                     Wrng=c(0,5), DWrng=c(-0.2, 0.2),
                     plotindex=1:n, 
                     titlestr = NULL, plotTitle = TRUE, xlab = "Score Index", ylab = NULL,
                     autoplot = FALSE,
                     plotMissing = TRUE, 
                     plotrange=c(min(scrfine),max(scrfine)), shaderange = NULL,  
                     ttlsz = NULL, axisttl = NULL, axistxt = NULL, 
                     lgdlab = NULL, lgdpos = "bottom") {
  # Plot a list of ICC plots
  
  # scrfine    Vector of length 101 containing plotting points
  # WfdList    List of Wbinsmth results
  # dataList   List of make.dataList
  # Qvec       vector of quantile positions
  
  # binctr     vector of bin centers, for plotting data points and/or confidence interval
  # data_point default as FALSE; set to TRUE for plotting data points
  # ci         default as FALSE; set to TRUE for plotting confidence limits of probability and surprisal
  
  # plotType   Type(s) of plot, default as "P" for probability, can also be "W" for surprisal, 
  #            "DW" for sensitivity, and any combination of the three
  
  # Wrng       make ylim of all suprisal plots the same
  # DWrng      make ylim of all sensitivity plots the same
  
  # plotindex  a vector of item indices; set if uses want to focus on specific items 
  
  # titlestr   plot title
  # plotTitle  indicator of showing the plot title, default as TRUE
  # xlab       title of x-axis, default as "Score Index" 
  # ylab       title of y-axis, default as "Probability", Note: should work with plotType to define ylab for each subplot
  
  # autoplot   plot all items in a batch
  
  # plotMissing Determine if plot the extra option for missing/spoiled responses.
  
  # plotrange  a vector of length 2; set if users want to focus on specific score range 
  # shaderange a list of length 2 vector(s); set if users want to gray out specific score range(s)
  
  # ttlsz      font size of plot title
  # axisttl    font size of axis titles
  # axistxt    font size of axis labels
  # lgdlab     font size of legend labels
  # lgdpos     legend position, could be set as "None" to remove the legend
  
  # Last modified July 10 2023 by Juan Li
  
  n <- length(WfdList)
  
  if (!is.null(shaderange)) {
    if (is.list(shaderange))
    {
      shadeValid <- rep(TRUE, length(shaderange))
      for (ishade in 1:length(shaderange))
      {
        rng <- shaderange[[ishade]]
        if (rng[1] >= rng[2])
        {
          shadeValid[ishade] <- FALSE
        } else if (rng[1] < plotrange[1])
        {
          shaderange[[ishade]][1] = plotrange[1]
        } else if (rng[2] > plotrange[2])
        {
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

  #   Invoke plotCore that actually controls the plotting
  plotlist <- list() # 2022-06-17
  if (autoplot | length(plotindex) == 1) { # 2023-05-28
    #  plot all items in a batch
    for (i in seq_len(length(plotindex)))# 2023-05-28
    {
      iplot <- plotindex[i]
      p <- plotCore(iplot, scrfine, WfdList, dataList, Qvec, 
                    binctr, data_point, ci, 
                    plotType,
                    Wrng, DWrng, 
                    titlestr, plotTitle, xlab, ylab,
                    plotMissing, 
                    plotrange, shaderangeV,  
                    ttlsz, axisttl, axistxt,
                    lgdlab, lgdpos)    
      #print(p)
      plotlist[[i]] <- p
    }
  } else {
    #  Usual case of user control of plotting
    nplot <- length(plotindex)
    jplot <- 1
    iplot <- plotindex[jplot]
    while (jplot >= 1 & jplot <= nplot) {
      p <- plotCore(iplot, scrfine, WfdList, dataList, Qvec, 
                    binctr, data_point, ci, 
                    plotType,
                    Wrng, DWrng, 
                    titlestr, plotTitle, xlab, ylab,
                    plotMissing, 
                    plotrange, shaderangeV,  
                    ttlsz, axisttl, axistxt,
                    lgdlab, lgdpos)    
      #print(p)
      plotlist[[iplot]] <- p
      #  request next plot
      
      if (length(plotindex) > 1) {
        n1 <- readline(prompt=
                         paste("Press [N/n/enter] to next item;", 
                               "[P/p] to previous item;", 
                               "or position in plotindex: "))
        if (n1 == "N" | n1 == "n" | n1 == "")
        {
          jplot <- jplot + 1
          iplot <- plotindex[jplot]
        } else if (n1 == "P" | n1 == "p") {
          jplot <- jplot - 1
          iplot <- plotindex[jplot]
        } else {
          n1 <- suppressWarnings(as.integer(n1))
          if (!is.na(n1) & (n1 >= 1 & n1 <= nplot)){
            jplot <- n1
            iplot <- plotindex[jplot]
          } else {
            break
          }
        }
      } else {
        break
      }
      gc() # to  free up memory from the previous run    } # end iplot
    }
  } # end iplot
  
  return(plotlist)
}

