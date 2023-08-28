Wbinsmth_nom <- function(bdry_nom, WfdList_nom) {
  #  approximate probability, surprisal and 2 surprisal derivatives
  #  over 101 equally spaced values between values in bdry
  #  for the nominal model estimated by mirt.
  n       <- length(WfdList_nom)
  nfine   <- 101
  indfine <- seq(bdry_nom[1],bdry_nom[2],len=nfine)
  #  set up seven bspline basis functions of order 5.
  Wnbasis <- 7  #  number of basis functions
  Wnorder <- 5  #  Order of the basis functions
  Wbasis  <- create.bspline.basis(bdry_nom, Wnbasis, Wnorder)
  #  loop through Lists in WfdList_nom
  for (item in 1:n) {
    WListi   <- WfdList_nom[[item]]
    #  argument WfdList_nom already loaded with number of options M
    #  and nominal parameter matrix set up by mirt
    Mi      <- WListi$M
    # mirt parameter matrix for this item
    parmati <- WListi$Wfd
    #  compute exponentials of Mi linear basis functions
    expmati <- matrix(0,nfine,Mi)
    for (m in 1:Mi) expmati [,m] <- exp(parmati[m,1]*indfine + parmati[m,2])
    #  probability matrix
    probmati <- expmati/apply(expmati,1,sum)
    #  surprisal matrix
    surpmati <- -log(probmati)/log(Mi)
    #  approximate surprisal curves using spline basis
    Result <- smooth.basis(indfine, surpmati, Wbasis)
    Sfd <- Result$fd
    #  evaluate the first and second derivative matrices
    Dsurpmati  <- eval.fd(indfine, Sfd, 1)
    D2surpmati <- eval.fd(indfine, Sfd, 2)
    #  load these four matrices into WfdList_nom
    WListi$Pmatfine   <-   probmati
    WListi$Wmatfine   <-   surpmati
    WListi$DWmatfine  <-  Dsurpmati
    WListi$D2Wmatfine <- D2surpmati
    #  compute arclength mesh and length
    infovec           <- pracma::cumtrapz(indfine, sqrt(apply(Dsurpmati^2,1,sum)))
    WListi$infovec    <- infovec
    WListi$arclength  <- max(infovec)
    #  load the struct WList into List array
    WfdList_nom[[item]] <- WListi
  }
  return(WfdList_nom)
}
