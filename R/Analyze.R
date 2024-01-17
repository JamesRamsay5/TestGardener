Analyze <- function(index, indexQnt, dataList, ncycle = 10, itdisp = FALSE, 
                    verbose = FALSE) {
  
  #' This function analyses a set of data by cycling ncycle numbers of 
  #' times between estimating probability and surprisal curves and 
  #' finding the optimum score index value for each person.
  #' The information about the data is stored in the list object dataStr.
  
  #. Arguments:
  #. index    ... A vector of length N that contains the initial estimates
  #.              of the score index values for the persons.  Usually, 
  #               the initial vector is a set of N equally spaced values
  #               spanning the interval [0,100].
  #. indexQnt ... The bin boundaries separated by the bin centers over [0,100].
  #.              The boundaries are chosen so that the numbers of persons
  #               in the bins are roughly equal.
  #. dataList ... A named list object containing all the information about
  #               the data that is required for analysis and subsequent 
  #               displays and tables.
  #. ncycle   ... The number of cycles in the analysis.
  #. itdisp   ... A logical value that determines whether the sequence of
  #               iteratiOns in the person scores is display in each cycle.
  #. verbose  ... A logical value that determines whether severalresults  
  #               within each cycle are displayed. 

  # Last modified 17 January 2024 by Jim Ramsay

  nbin    <- dataList$nbin            # number of bins
  markers <- dataList$PcntMarkers/100 # marker probabilities
  bdry0   <- seq(0,2*nbin,1)/(2*nbin) # initial boundary values
  
  #  set up list vector to contain all results for each cycle
  
  parmListvec <- vector("list",ncycle)  
  
  #  define the spline basis for representing the log density function
  
  logdensbasis <- fda::create.bspline.basis(c(0,100), 15)
  
  #  extract information about surprisal smoothing for each item
  
  SfdList <- dataList$SfdList
  n       <- length(SfdList)
  
  #  compute the dimension of the space within which the surprisal curves
  #  vary
  
  Sdim <- 0
  for (i in 1:n) {
    SStri <- SfdList[[i]]
    Sdim  <- Sdim + SStri$M
  }
  
  #  extract the data matrix from argument dataList
  
  chcemat <- dataList$chcemat
  
  #. define a mesh of 101 score index values spanning interval [0,100]
  
  indfine   <- seq(0,100,len=101)
  
  #  ----------------------------------------------------------
  #                       main cycle loop
  #  ----------------------------------------------------------
  
  for (icycle in 1:ncycle) {
    
    print(paste('----------  Cycle ',icycle,'-----------'))
    
    #  ----------------------------------------------------------
    #  Step 1:  Bin the data, and smooth the binned data
    #  ----------------------------------------------------------
    
    # print("step 1")
    
    if (verbose) print("Optimize surprisal curves:")
    
    #  Sbinsmth uses bin boundaries and centres in argument indexQnt
    #  to allocate score indices to bins, compute proportions their
    #  surprisal values, and then loop through items to estimate
    #  smooth probability and surprisal curves.
    #  Function smooth.surp is used for this purpose.
    #  After score index values are computed, bin boundaries and 
    #  and centres are adjusted in Step 4.
    
    SfdResult <- Sbinsmth(index, dataList)
    SfdList   <- SfdResult$SfdList
    binctr    <- SfdResult$binctr
    bdry      <- SfdResult$bdry
    freq      <- SfdResult$freq
    
    # print("Step 1 finished")
    
    # print("bin frequencies:")
    # print(t(freq))
    # print("bin ctrs:")
    # print(round(binctr,1))
          
    #  ----------------------------------------------------------
    #  Step 2:  compute mean value of objective function 
    #  ----------------------------------------------------------
    
    # print("step 2")
    
    if (verbose) print("Compute mean examinee fits")
    
    Fvec  <- Ffun(index, SfdList, chcemat)
    meanF <- mean(Fvec)
    
    if (verbose) print(paste('Mean data fit = ', round(meanF,3)))
    
    # print("Step 2 finished")
    
    #  ------------------------------------------------------------
    #  Step 3:  Compute optimal score index values for each person
    #  ------------------------------------------------------------
    
    # print("step 3")
    
    # if (verbose) print("Optimize examinee data fits")

    # optimize score index values
    
    indexfunList <- index_fun(index, SfdList, chcemat, 20, 1e-3, itdisp=itdisp)
    
    # extract information
    
    index    <- indexfunList$index_out
    Fval     <- indexfunList$Fval
    DFval    <- indexfunList$DFval
    D2Fval   <- indexfunList$D2Fval
    active   <- indexfunList$active
    
    # print("Step 3 finished")
    
    #  ----------------------------------------------------------
    #  Step 4:  Estimate the score density for score index values
    #  The density is only defined by score index values inside
    #  [0,100].  Counts of values on the boundaries are indicated
    #  circles on the boundary.  The density is used to adjust
    #  bin boundaries and centres for the next cycle.
    #  Step 4 is not taken if this is the final cycle.
    #  ----------------------------------------------------------
    
    if (icycle < ncycle) {
    # print("step 4")
    
    # if (verbose) print("Compute score index density")
    
    # use only interior score index values
    indexdens <- index[0 < index & index < 100]
    # estimate the cumulative density denscdf over values indcdf
    index_distnList <- index_distn(indexdens, logdensbasis)
    denscdf         <- as.numeric(index_distnList$denscdf)
    indcdf          <- as.numeric(index_distnList$indcdf)
    
    # # adjusted marker score index values are computed by interpolation
    Qvec <- pracma::interp1(denscdf, indcdf, markers)
    
    density_plot(indexdens, c(0,100), Qvec, xlabstr="Score index",
                 titlestr=paste("Current index density, cycle",icycle),
                 scrnbasis=11, nfine=101)
    
    # This interpolation adjusts bin boundaries and centres to define
    # a new vector indexQnt
    # print(round(denscdf,3))
    # print(round(indcdf,3))
    
    #. compute 2*nbin - 1 inner boundary/center pair locations by interpolation
    
    # innerindex <- seq(2,2*nbin,1)
    # innerQnt   <- indexQnt[innerindex]
    # innerQnt   <- pracma::interp1(as.numeric(denscdf), as.numeric(indcdf/100),
    #                             innerQnt/100)*100
    # #. define new version of complete indexQnt
    # indexQnt[innerindex] <- innerQnt
    
    #. bin centres
    # plot(indcdf, denscdf, type="b", xlim=c(0,100), ylim=c(0,1))
    # for (i in seq(1,2*nbin+1,2)) lines(c(indexQnt[i],indexQnt[i]), c(0,1))
    binctr <- indexQnt[seq(2,2*nbin,  2)]
    #. bin boundaries
    bdry   <- indexQnt[seq(1,2*nbin+1,2)]
    # print("Bin boundaries")
    # print(round(bdry,1))
    # print("Bin centres")
    # print(round(binctr,1))

    # readline(prompt = "Enter to continue:")

    # print("Step 4 finished")
      
    }
    
    #  ----------------------------------------------------------
    #  Step 5.  Compute arc length of the surprisal space curve
    #  ----------------------------------------------------------
    
    # print("step 5")
    
    DSfine <- matrix(0,101,Sdim)
    m2 <- 0
    for (i in 1:n) {
      SListi <- SfdList[[i]]
      Mi     <- SListi$M
      m1 <- m2 + 1
      m2 <- m2 + Mi
      DSfine[,m1:m2] <- SListi$DSmatfine
    }
    indfine  <- seq(0,100,len=101)
    infoSurp <- max(pracma::cumtrapz(indfine, sqrt(apply(DSfine^2,1,sum))))
    
    if (verbose)  print(paste('infoSurp in bits = ',round(infoSurp,1)))
    
    # print("Step 5 finished")
    
    #  ----------------------------------------------------------
    #  Step 6:  Check for mis-identifications of minimum index 
    #  ----------------------------------------------------------
    
    # print("step 6")
    
    Result <- index_search(SfdList, dataList$chcemat, index, 
                             Fval, DFval, D2Fval)
    index  <- Result$index
    Fval   <- Result$Fval
    DFval  <- Result$DFval
    D2Fval <- Result$D2Fval
    
    # print("Step 6 finished")
    
    #  -----------------------------------------------------------
    #  Step 7:  set up the parameter list parmListi for this cycle
    #  -----------------------------------------------------------

    # print("step 7")
    
    parmListi <- list(
      index      = index,
      indexQnt   = indexQnt,
      SfdList    = SfdList,
      meanF      = meanF,
      binctr     = binctr,
      bdry       = bdry,
      freq       = freq,
      Qvec       = Qvec,
      Fval       = Fval,
      DFval      = DFval,
      D2Fval     = D2Fval,
      active     = active,
      infoSurp   = infoSurp
    )
    
    parmListvec[[icycle]] <- parmListi
    
    # print("Step 7 finished")
    
  }
  
  # end of the loop through the cycles.  
  
  #. return the list object parmListvec of length ncycle
  #. containing parameter results for each cycle
  
  return(parmListvec=parmListvec)
  
}