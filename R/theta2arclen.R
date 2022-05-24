theta2arclen <- function(theta, Qvec, WfdList, itemindex=1:n, plotrng=c(0,100)) {
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
    # WFDLIST   ... Cell array containing estimated surprisal info.
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
    # WFD.INFO      ... log derivative functional data object defining
    #                   a strictly increasing set of score index values
    #                   corresponding to set of arc length values
    # THETAVEC      ... score index values resulting from using function
    #                   monfd with equally spaced arc length values
    #                   and WFD.INFO.
    # WDIM          ... The dimension of the over space containing the 
    #                   surprisal curves.
    
    #  Last modified 16 May 2022
    
    #  ------------------------------------------------------------------------
    #                              check inputs  
    #  ------------------------------------------------------------------------
    
    if (!is.list(WfdList)) {
        stop("The second argument is not a cell object.")
    }
    
    if (!is.numeric(theta)) {
        stop("The first argument is not of class numeric.")
    }
    
    if (any(theta < 0) || any(theta > 100)) {
        stop("The first argument has values outside of c(0,100).")
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
    
    indfine <- seq(plotrng[1], plotrng[2], len=101)
    nfine   <- length(indfine)
    
    #  compute Wdim, the dimension of the overspace within which the
    #  surprisal curves vary
    
    Wdim <- 0
    for (i in itemindex) {
        WListi <- WfdList[[i]]
        Wdim  <- Wdim + WListi$M
    }
    
    #  set up surprisal derivatives over a fine mesh
    
    DWfine <- matrix(0,nfine,Wdim)
    m2 <- 0
    for (i in itemindex) {
        WListi <- WfdList[[i]]
        Mi     <- WListi$M
        m1 <- m2 + 1
        m2 <- m2 + Mi
        DWfine[,m1:m2] <- WListi$DWmatfine
    }
    
    #  ------------------------------------------------------------------------
    #  Use the trapezoidal rule to define a sequence of length 1001
    #  over the surprisal manifold defined by surprisal curve values 
    #  associated with equally spaced values over the score index continuum.
    #  This sequence us not equally spaced over the arc length continuum,
    #  however.
    #  Trapezoidal rule indefinite integration over 101 points is not 
    #  very accurate, but given the smoothness of the relation between
    #  theta and arclength or information, this should be adequate for 
    #  plotting purposes.
    #  A more accurate method using Richard extrapolation might be needed
    #  at some point.
    #  ------------------------------------------------------------------------
    
    arclengthvec <- pracma::cumtrapz(sqrt(apply(DWfine^2,1,sum)))
    arclength    <- max(arclengthvec)
    
    #  ------------------------------------------------------------------------
    #  Compute arc length values for equally spaced theta values and 
    #  marker percentages using interpolation over  
    #  the fine mesh of score index values.
    #  ------------------------------------------------------------------------
    
    Wnbasis <- 11
    
    Wbasis.theta <- create.bspline.basis(c(0,100),Wnbasis)
    WfdPar.theta <- fdPar(fd(matrix(0,Wnbasis,1),Wbasis.theta))
    Wfd.theta    <- smooth.morph(indfine, arclengthvec, c(0,arclength), WfdPar.theta)$Wfd
    monfnmax     <- monfn(100, Wfd.theta)
    
    #  compute arc length values corresponding to N estimated theta values
    
    theta_al <- monfn(theta, Wfd.theta)*arclength/monfnmax
    
    #  compute arc length values corresponding to estimated marker theta values
    
    Qvec_al  <- monfn(Qvec, Wfd.theta)*arclength/monfnmax
    
    #  ------------------------------------------------------------------------
    # Also compute 101 theta values that yield a roughly equally spaced set of
    # arc length values. 
    #  ------------------------------------------------------------------------
    
    Wbasis.info <- create.bspline.basis(c(0,arclength), Wnbasis)
    WfdPar.info <- fdPar(fd(matrix(0,Wnbasis,1),Wbasis.info))
    Wfd.info    <- smooth.morph(arclengthvec, indfine, c(0,100), WfdPar.info)$Wfd
    monfnmax    <- monfn(arclength, Wfd.info)
    thetavec    <- monfn(arclengthvec, Wfd.info)*100/monfnmax
    
    #  return a list object infoList to contain all the objects
    
    return(list(
        arclength     = arclength,
        Wfd.theta     = Wfd.theta,
        arclengthvec  = arclengthvec,
        theta_al      = theta_al,
        Qvec_al       = Qvec_al,
        Wfd.info      = Wfd.info,
        thetavec      = thetavec,
        Wdim          = Wdim))
    
}