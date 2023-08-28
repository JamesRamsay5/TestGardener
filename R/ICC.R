
ICC <- function(x) {
  UseMethod("ICC")
}

ICC <- function(Item, M, Sfd, Pbin, Sbin, Pmatfine, Sarrayfine, 
                PStdErr, SStdErr, ItemArcLen, itemStr=NULL, optStr=NULL) {
  # ICC is the container of information required for an item and 
  # the options within it.
  # Arguments:
  # Item.      ...  The index of the item within the sequence 1,...,n
  # M          ...  Length of probability and surprisal vector
  # Sfd        ...  A functional data object whose curves satisfy the 
  #                 surprisal constraints.
  # Pbin       ...  An nbin by M matrix of probabilities, summing to one 
  #                 across columns. The number of rows is the number of 
  #                 bins for estimating probabilities of choice, and M is 
  #                 the number of options.
  # Sbin       ...  An nbin by M matrix of surpisal values, 
  #                 surprisal = -log_M(probability).
  # Pmatfine   ...  A 101 by M probability values over an equally spaced 
  #                 mesh length 101.
  # Sarrayfine ...  A 101 by M by 3 array containing suprisal values over 
  #                 an equally spaced mesh of size 101 in the first layer,
  #                 surprisal derivatives in the second and  surprisal
  #                 second drivatives in the third
  # PStdErr    ...  Standard errors at bin centers for probability curves
  # SStdErr    ...  Standard errors at bin centers for surprisal curves
  # ItemArcLen ...  Arc length for item information curve
  # itemStr    ...  Item label string
  # optStr     ...  List vector containing option label strings
  
  #  Last modified 20 April 2023 by Jim Ramsay
  
  #  Set up default binary ICC if there are no arguments
  
  if (nargs()==0) { 
    Item   <-  1
    M      <-  2
    nbin   <- 11
    nbasis <-  3
    norder <-  3
    basisobj   <- create.bspline.basis(c(0,100), nbasis, norder)
    Sfd        <- fd(matrix(1,nbasis,M), basisobj)
    Pbin       <- matrix(1/M,nbin,M)
    Sbin       <- matrix(  1,nbin,M)
    Pmatfine   <- matrix(1/M, 101,M)
    Sarrayfine <- array(0, c(101,M,3))
    PStdErr    <- matrix(0,nbin)
    SStdErr    <- matrix(0,nbin)
    ItemArcLen <- 0
    itemStr    <- NULL
    optStr     <- NULL
    ICCList <- list(Item=Item,         Sfd=Sfd, 
                    Pbin=Pbin,         Sbin=Sbin, 
                    Pmatfine=Pmatfine, Sarrayfine=Sarrayfine,
                    PStdErr=PStdErr,   PStdErr=PStdErr,
                    ItemArcLen=ItemArcLen, 
                    itemStr=itemStr,  optStr=optStr)
    class(ICC) <- "ICC"
    return(ICC)
  }
  
  #  if first argument is a basis object, return
  
  if (inherits(ICC, "ICC")) {
    return(ICC)
  }
  
  # check Item
  
  if (!is.numeric(Item)) 
    stop("Item is not a numeric object.  It must be the index of the item.")
  
  # check M
  
  if (!is.integer(M)) stop("M is not an integer.")
  if (!(M > 1)) stop("M is not greater than 1.")
  
  # check if Sfd is an fd object and has a bspline basis
  
  if (inherits(Sfd, "fd")) {
    if (ncol(Sfd$coefs) != M) stop("ncol(Sfd$coefs) is not M.")
    if (!inherits(Sfd$basis, "basisfd")) 
      stop("Sfd$basis is not a basis object.")
    if (!Sfd$basis$type == "bspline")
      stop("Sfd does not have a bspline basis.")
  } else {
    stop("Sfd is not a functional data object.")
  }
  
  # check Pbin and Sbin
  
  if (is.matrix(Pbin) && is.matrix(Sbin)) {
    if (ncol(Pbin) == M && ncol(Sbin) == M) {
      nbin <- nrow(Pbin)
      if (nbin != nrow(Sbin)) 
        stop("Pbin and Sbin have different numbers of rows.")
    } else {
      stop(paste("Number of columns of Pbin and/or Sbin not equal",
                 "to ncol(Sfd$coefs)."))
    }
  } else {
    stop("Pbin and/or Sbin are not matrices")
  }
  
  # check Pmatfine and Sarrayfine
  
  if (!is.matrix(Pmatfine)) stop("Pmatfine is not a matrix.")
  if (!(ncol(Pmatfine) == M && nrow(Pmatfine) == 101)) 
        stop("Pmatfine does not have M columns and 101 rows.")
  
  if (!is.array(Sarrayfine)) stop("Sarrayfine is not an array.")
  Sdim <- dim(Sarrayfine)
  if (!(Sdim[1] == 101 && Sdim[2] == M && Sdim[3] == 3)) 
    stop("Sarrayfine does not have dimensions 101, M and 3.")

  #  check PStdErr and SStdErr
  
  if (!is.matrix(PStdErr)) stop("PStdErr is not a matrix.")
  if (!(ncol(PStdErr) == M && nrow(PStdErr) == nbin)) 
    stop("PStdErr does not have M columns and nbin rows.")
  
  if (!is.matrix(SStdErr)) stop("SStdErr is not a matrix.")
  if (!(ncol(SStdErr) == M && nrow(SStdErr) == nbin)) 
    stop("SStdErr does not have M columns and nbin rows.")
  
  #  check ItemArcLen
  
  if (!is.numeric(ItemArcLen)) stop("ItemArcLen is not numeric.")
  if (ItemArcLen < 0) stop("ItemArcLen is negative.")
  
  ICCList <- list(Item=Item,        Sfd=Sfd, 
                  Pbin=Pbin,         Sbin=Sbin, 
                  Pmatfine=Pmatfine, Sarrayfine=Sarrayfine,
                  PStdErr=PStdErr,   PStdErr=PStdErr,
                  ItemArcLen=ItemArcLen, 
                  itemStr=itemStr,  optStr=optStr)
  
  class(ICC) <- "ICC"
  
  return(ICCList)
  
}

#  ----------------------------------------------------------------------------

print.ICC <- function(x, Pbin, Sbin, ...) {
  print("Binned Probabilities:")
  print(round(Pbin,3))
  print("Binned Surprisals:")
  print(round(Sbin,3))
}

#  ----------------------------------------------------------------------------

summary.ICC <- function(x, ...) {
  print("Sfd: Surprisal functional data object for surprisal curves.")
  print("Pbin: nbin by M matrix of proportions of choices.")
  print("Sbin: nbin by M matrix of surprisal values of choices.")
  print("Pmatfine: 101 by M matrix of probability curve values.")
  print(paste("Sarray: 101 by M by 3 array of suprisal values", 
              "and their first 2 derivatives"))
  print(paste("itemStr: item label string"))
  print(paste("optStr:  list vector of option strings."))
}

  