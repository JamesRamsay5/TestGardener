index2info <- function(index, Qvec, SfdList, binctr, itemindex=1:n, 
                      plotrng=c(0,100), shortwrd=FALSE) {
    # INDEX2INFO is an important function.  It's job is to convert
    #  objects defined over the score index continuum [0,100] to
    # the same objects over the arc length continuum [0,infoSurp], and also
    # vice versa.  Since the arc length or information continuum is along
    # a space curve that is invariant under strictly monotone transformations
    # of the score index \index, and is also a metric, it is an ideal 
    #  choice for the abscissa in all plots.
    #
    #  It is called in functions:
    #    
    # TestGardener/analysiscode/Analyze.m
    #
    # Arguments:
    #
    # INDEX     ... A single column vector of score index values for each of
    #               N test takers.  The values are over {0,100].
    # QVEC      ... Estimated marker percentages for score index.
    # SFDLIST   ... List array containing estimated surprisal info.
    # BINCTR    ... Bin centers for score index
    # ITEMINDEX ... Vector of item indices to use.
    # PLOTRNG   ... Two-vector specifying the score index range to use.
    #
    # Returned objects:
    #
    # INFOSURP      ... arc length = max(int.infoSurp)
    # SFD.INDEX     ... log derivative functional data object defining
    #                   a strictly increasing set of arc length values
    #                   corresponding to set of score index values
    # INFOSURPVEC   ... values of indefinite integral of sum of norms
    #                   of surprisal derivatives using trapezoidal rule
    # SCOPEVEC      ... arc length values corresponding to estimated 
    #                   score index value
    # QINFOVEC      ... arc length marker percentages corresponding to
    #                   marker percentages for score index values
    # BININFOCTR    ... Bin centers for arc length or information
    # SFD.INFO      ... log derivative functional data object defining
    #                   a strictly increasing set of score index values
    #                   corresponding to set of arc length values
    # SDIM          ... The dimension of the over space containing the 
    #                   surprisal curves.
    
    #  Last modified 5 January 2024
    
    #  ------------------------------------------------------------------------
    #                              check inputs  
    #  ------------------------------------------------------------------------
    
    if (!is.list(SfdList)) {
        stop("The third argument is not a list object.")
    }
    
    if (!is.numeric(index)) {
        stop("The first argument is not of class numeric.")
    }
  
    if (any(index < plotrng[1]) || any(index > plotrng[2])) {
        stop("The first argument has values outside of argument plotrng.")
    }
    
    #  ------------------------------------------------------------------------
    #  define number of items n, assign default values to last two 
    #  arguments and check the last two arguments
    #  ------------------------------------------------------------------------
    
    n <- length(SfdList)
    
    if (!is.numeric(plotrng) || !is.numeric(itemindex)) {
        stop("Arguments plotrng or itemindex are not numeric.")
    } 
    
    if (is.null(itemindex)) {
        stop("Argument itemindex is empty.")
    }
    
    if (length(plotrng) != 2) {
        stop("Argument plotrng is not of length 2.")
    }
    
    if (min(itemindex) < 1 || max(itemindex) > n+1) {
        stop("Values in itemindex are < 1 or > n+1.")
    }
    
    #  ------------------------------------------------------------------------
    #  Set up a fine mesh over score index range defined in argument plotrng
    #  Here and elsewhere 'fine' implies equal spacing.
    #  ------------------------------------------------------------------------
    
    index <- index[index >= plotrng[1] & index <= plotrng[2]]
    nfine <- 101
    indfine.rng <- as.matrix(seq(plotrng[1], plotrng[2], len=nfine))
    
    #  ------------------------------------------------------------------------
    #  compute Sdim, the dimension of the overspace within which the
    #  surprisal curves vary
    #  ------------------------------------------------------------------------
    
    Sdim <- 0
    for (i in itemindex) {
        Sdim   <- Sdim + SfdList[[i]]$M
    }
    
    #  ------------------------------------------------------------------------
    #  set up surprisal derivatives over plot range
    #  ------------------------------------------------------------------------
    
    DSfine.rng <- matrix(0,nfine,Sdim)
    m2 <- 0
    for (item in itemindex) {
        SListi <- SfdList[[item]]
        Mi     <- SListi$M
        m1 <- m2 + 1
        m2 <- m2 + Mi
        DSfine.rng[,m1:m2] <- SListi$DSmatfine
    }

    #  ------------------------------------------------------------------------
    #  Compute arc length values infofine for equally spaced index values over  
    #  the items in plotindex over the plot range.  Integration is by the 
    #. trapezoidal rule.
    #  ------------------------------------------------------------------------
    
    infoSurpvec.rng <- 
      pracma::cumtrapz(indfine.rng,sqrt(apply(DSfine.rng^2,1,sum)))
    infoSurp.rng    <- max(infoSurpvec.rng)

    #  ------------------------------------------------------------------------
    #  If shortwrd == TRUE, return here with only infoSurp.rng and
    # infoSurpvec.rng
    #  ------------------------------------------------------------------------

    if (shortwrd) {
      return(list(infoSurp = infoSurp.rng,  infoSurpvec = infoSurpvec.rng))
    }   
    
    #  ------------------------------------------------------------------------
    #  set up a fd object over plot range for representing arc length
    #  ------------------------------------------------------------------------
    
    Snbasis.rng <- 11
    Sbasis.rng  <- fda::create.bspline.basis(plotrng,Snbasis.rng)
    Sfd.rng     <- fda::fd(matrix(0,Snbasis.rng,1), Sbasis.rng)
    
    #  ------------------------------------------------------------------------
    #  Compute the monotone functional data object transforming the score index 
    #  values to arc length values over the plotting range.
    #  ------------------------------------------------------------------------
    
    infoSurpfd.rng <- fda::smooth.morph(indfine.rng, infoSurpvec.rng, 
                                   c(0,infoSurp.rng), Sfd.rng)$Wfd
    #  evaluate the monotone function at arc length
    monfnmax        <- fda::monfn(infoSurp.rng, infoSurpfd.rng)
    
    #  ------------------------------------------------------------------------
    #  compute arc length values corresponding to estimated marker index values
    #  ------------------------------------------------------------------------
    
    Qvec.rng <- Qvec[Qvec >= plotrng[1] & Qvec <= plotrng[2]]
    if (length(Qvec.rng) > 0) {
      Qinfovec  <- pracma::interp1(as.numeric(indfine.rng), 
                                  as.numeric(infoSurpvec.rng), Qvec.rng)
    } else {
      Qinfovec <- NULL
    }
    
    #  ------------------------------------------------------------------------
    #  compute arc length values corresponding to bin centers
    #  ------------------------------------------------------------------------
    
    if (!is.null(binctr)) {
    binctr.rng <- binctr[binctr >= plotrng[1] & binctr <= plotrng[2]]
    bininfoctr  <- pracma::interp1(as.numeric(indfine.rng), 
                                as.numeric(infoSurpvec.rng), binctr.rng)
    } else {
      binctr.rng <- NULL
      bininfoctr  <- NULL
    }
    
    #  ------------------------------------------------------------------------
    #  Compute arc length values scopevec corresponding to estimated index 
    #  values within plot range
    #  ------------------------------------------------------------------------
    
    if (is.null(index)) {
      scopevec <- NULL
    } else {
      index.rng     <- index[index >= plotrng[1] & index <= plotrng[2]]
      scopevec.rng  <- pracma::interp1(as.numeric(indfine.rng), 
                                   as.numeric(infoSurpvec.rng), index.rng)
    }
    
    #  ------------------------------------------------------------------------
    # Also compute 101 index values that yield a roughly equally spaced set of
    # arc length values. 
    #  ------------------------------------------------------------------------
    
    Sbasis.info <- fda::create.bspline.basis(c(0,infoSurp.rng), Snbasis.rng)
    Sfd.info    <- fda::fd(matrix(0,Snbasis.rng,1),Sbasis.info)
    Sfd.info    <- fda::smooth.morph(infoSurpvec.rng, indfine.rng, plotrng, 
                                     Sfd.info)$Wfd
    monfnmax          <- fda::monfn(infoSurp.rng, Sfd.info)
    infoSurpfine.rng  <- seq(0,infoSurp.rng,leng=nfine)
    plotwidth         <- plotrng[2] - plotrng[1]
    infofine.rng      <- plotrng[1] + 
      fda::monfn(infoSurpfine.rng, Sfd.info)*plotwidth/monfnmax
    
    #  ------------------------------------------------------------------------
    #  return a list object infoList to contain all the objects
    #  ------------------------------------------------------------------------
    
    return(list(
        infoSurp    = infoSurp.rng,
        infoSurpvec = infoSurpvec.rng,
        infoSurpfd  = infoSurpfd.rng,
        scopevec    = scopevec.rng,
        Qinfovec    = Qinfovec,
        bininfoctr  = bininfoctr,
        Sfd.info    = Sfd.info,
        Sdim        = Sdim))
    
}