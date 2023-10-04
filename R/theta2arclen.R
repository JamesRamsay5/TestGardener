theta2arclen <- function(theta, Qvec, WfdList, binctr, itemindex=1:n, 
                         plotrng=c(0,100), shortwrd=FALSE) {
    # THETA2ARCLEN is a centrally important function.  It's job is to convert
    #  objects defined over the score index continuum [0,100] to
    # the same objects over the arc length continuum [0,arclength], and also
    # vice versa.  Since the arc length or information continuum is along
    # a space curve that is invariant under strictly monotone transformations
    # of the score index \theta, and is also a metric, it is an ideal 
    #  choice for the abscissa in all plots.
    #
    #  It is called in functions:
    #    
    # TestGardener/analysiscode/analyze.m
    # TestGardener/analysiscode/Conditional.simulaton.m
    # TestGardener/analysiscode/Inforng.plot
    #
    # Arguments:
    #
    # THETA     ... A single column vector of score index values for each of
    #               N test takers.  The values are over {0,100].
    # QVEC      ... Estimated marker percentages for score index.
    # WFDLIST   ... List array containing estimated surprisal info.
    # BINCTR    ... Bin centers for score index
    # ITEMINDEX ... Vector of item indices to use.
    # PLOTRNG   ... Two-vector specifying the score index range to use.
    #
    # Returned objects:
    #
    # ARCLENGTH     ... arc length = max(int.arclength)
    # WFD.THETA     ... log derivative functional data object defining
    #                   a strictly increasing set of arc length values
    #                   corresponding to set of score index values
    # ARCLENGTHVEC  ... values of indefinite integral of sum of norms
    #                   of surprisal derivatives using trapezoidal rule
    # THETA.AL      ... arc length values corresponding to estimated 
    #                   score index value
    # QVEC.AL       ... arc length marker percentages corresponding to
    #                   marker percentages for score index values
    # BINCTR.AL     ... Bin centers for arc length or information
    # WFD.INFO      ... log derivative functional data object defining
    #                   a strictly increasing set of score index values
    #                   corresponding to set of arc length values
    # THETAVEC      ... score index values resulting from using function
    #                   monfd with equally spaced arc length values
    #                   and WFD.INFO.
    # WDIM          ... The dimension of the over space containing the 
    #                   surprisal curves.
    
    #  Last modified 25 September 2023
    
    #  ------------------------------------------------------------------------
    #                              check inputs  
    #  ------------------------------------------------------------------------
    
    if (!is.list(WfdList)) {
        stop("The third argument is not a list object.")
    }
    
    if (!is.numeric(theta)) {
        stop("The first argument is not of class numeric.")
    }
  
    if (any(theta < plotrng[1]) || any(theta > plotrng[2])) {
        stop("The first argument has values outside of argument plotrng.")
    }
    
    #  ------------------------------------------------------------------------
    #  define number of items n, assign default values to last two 
    #  arguments and check the last two arguments
    #  ------------------------------------------------------------------------
    
    n <- length(WfdList)
    
    if (!is.numeric(plotrng) || !is.numeric(itemindex)) {
        stop("Arguments plotrng or itemindex are not numeric.")
    } 
    
    if (is.null(itemindex)) {
        stop("Argument itemindex is empty.")
    }
    
    if (length(plotrng) != 2) {
        stop("Argument plotrng is not of length 2.")
    }
    
    if (min(itemindex) < 1 || max(itemindex) > n) {
        stop("Values in itemindex are < 1 or > n.")
    }
    
    #  ------------------------------------------------------------------------
    #  Set up a fine mesh over score index range defined in argument plotrng
    #  Here and elsewhere 'fine' implies equal spacing.
    #  ------------------------------------------------------------------------
    
    theta <- theta[theta >= plotrng[1] & theta <= plotrng[2]]
    nfine <- 101
    indfine.rng <- as.matrix(seq(plotrng[1], plotrng[2], len=nfine))
    
    #  ------------------------------------------------------------------------
    #  compute Wdim, the dimension of the overspace within which the
    #  surprisal curves vary
    #  ------------------------------------------------------------------------
    
    Wdim <- 0
    for (i in itemindex) {
        Wdim   <- Wdim + WfdList[[i]]$M
    }
    
    #  ------------------------------------------------------------------------
    #  set up surprisal derivatives over plot range
    #  ------------------------------------------------------------------------
    
    DWfine.rng <- matrix(0,nfine,Wdim)
    m2 <- 0
    for (item in itemindex) {
        WListi <- WfdList[[item]]
        Mi     <- WListi$M
        m1 <- m2 + 1
        m2 <- m2 + Mi
        DWfine.rng[,m1:m2] <- WListi$DWmatfine
    }

    #  ------------------------------------------------------------------------
    #  Compute arc length values for equally spaced theta values over the 
    #  items in plotindex over the plot range.  Integration is by the 
    #. trapezoidal rule.
    #  ------------------------------------------------------------------------
    
    arclengthvec.rng <- 
      pracma::cumtrapz(indfine.rng,sqrt(apply(DWfine.rng^2,1,sum)))
    arclength.rng    <- max(arclengthvec.rng)

    scopevec.rng <- matrix(0,length(itemindex),1)
    for (item in itemindex) {
      WListi <- WfdList[[item]]
      scopevec.rngi <- 
        pracma::cumtrapz(indfine.rng,sqrt(apply(WListi$DWmatfine^2,1,sum)))
      scopevec.rng[item]    <- max(scopevec.rngi)
    }
    
    #  ------------------------------------------------------------------------
    #  If shortwrd == TRUE, return here with only arclength.rng and
    # arclengthvec.rng
    #  ------------------------------------------------------------------------

    if (shortwrd) {
      return(list(arclength = arclength.rng,  arclengthvec = arclengthvec.rng))
    }   
    
    #  ------------------------------------------------------------------------
    #  set up a fd object over plot range for representing arc length
    #  ------------------------------------------------------------------------
    
    Wnbasis.rng <- 11
    Wbasis.rng  <- create.bspline.basis(plotrng,Wnbasis.rng)
    Wfd.rng     <- fd(matrix(0,Wnbasis.rng,1), Wbasis.rng)
    WfdPar.rng  <- fdPar(Wfd.rng)
    
    #  ------------------------------------------------------------------------
    #  Compute the monotone functional data object transforming the score index 
    #  values to arc length values over the plotting range.
    #  ------------------------------------------------------------------------
    
    arclengthfd.rng <- smooth.morph(indfine.rng, arclengthvec.rng, 
                                    c(0,arclength.rng), WfdPar.rng)$Wfd
    #  evaluate the monotone function at arc length
    monfnmax        <- monfn(arclength.rng, arclengthfd.rng)
    
    #  ------------------------------------------------------------------------
    #  compute arc length values corresponding to estimated marker theta values
    #  ------------------------------------------------------------------------
    
    Qvec.rng <- Qvec[Qvec >= plotrng[1] & Qvec <= plotrng[2]]
    if (length(Qvec.rng) > 0) {
      Qvec_al  <- pracma::interp1(as.numeric(indfine.rng), 
                                  as.numeric(arclengthvec.rng), Qvec.rng)
    } else {
      Qvec_al <- NULL
    }
    
    #  ------------------------------------------------------------------------
    #  compute arc length values corresponding to bin centers
    #  ------------------------------------------------------------------------
    
    if (!is.null(binctr)) {
    binctr.rng <- binctr[binctr >= plotrng[1] & binctr <= plotrng[2]]
    binctr_al  <- pracma::interp1(as.numeric(indfine.rng), 
                                as.numeric(arclengthvec.rng), binctr.rng)
    } else {
      binctr.rng <- NULL
      binctr_al  <- NULL
    }
    
    #  ------------------------------------------------------------------------
    #  compute arc length values corresponding to estimated theta values
    #  within plot range
    #  ------------------------------------------------------------------------
    
    if (is.null(theta)) {
      theta_al <- NULL
    } else {
      theta.rng <- theta[theta >= plotrng[1] & theta <= plotrng[2]]
      theta_al  <- pracma::interp1(as.numeric(indfine.rng), 
                                    as.numeric(arclengthvec.rng), theta.rng)
    }
    
    #  ------------------------------------------------------------------------
    # Also compute 101 theta values that yield a roughly equally spaced set of
    # arc length values. 
    #  ------------------------------------------------------------------------
    
    Wbasis.info <- create.bspline.basis(c(0,arclength.rng), Wnbasis.rng)
    WfdPar.info <- fdPar(fd(matrix(0,Wnbasis.rng,1),Wbasis.info))
    Wfd.info    <- smooth.morph(arclengthvec.rng, indfine.rng, plotrng, 
                                WfdPar.info)$Wfd
    monfnmax          <- monfn(arclength.rng, Wfd.info)
    arclengthfine.rng <- seq(0,arclength.rng,leng=nfine)
    plotwidth         <- plotrng[2] - plotrng[1]
    thetafine.rng     <- plotrng[1] + 
                         monfn(arclengthfine.rng, Wfd.info)*plotwidth/monfnmax
    
    #  ------------------------------------------------------------------------
    #  return a list object infoList to contain all the objects
    #  ------------------------------------------------------------------------
    
    return(list(
        arclength       = arclength.rng,
        scopevec.       = scopevec.rng,
        arclengthvec    = arclengthvec.rng,
        arclengthfd     = arclengthfd.rng,
        theta_al        = theta_al,
        thetafine_al    = thetafine.rng,
        Qvec_al         = Qvec_al,
        binctr_al       = binctr_al,
        Wfd.info        = Wfd.info,
        Wdim            = Wdim))
    
}