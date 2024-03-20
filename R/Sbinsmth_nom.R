Sbinsmth_nom <- function(bdry_nom, SfdList_nom) {
  #  approximate probability, surprisal and 2 surprisal derivatives
  #  over 101 equally spaced values between values in bdry
  #  for the nominal model estimated by mirt.
  n       <- length(SfdList_nom)
  nfine   <- 101
  indfine <- seq(bdry_nom[1],bdry_nom[2],len=nfine)
  #  set up seven bspline basis functions of order 5.
  Snbasis <- 7  #  number of basis functions
  Snorder <- 5  #  Order of the basis functions
  Sbasis  <- create.bspline.basis(bdry_nom, Snbasis, Snorder)
  #  loop through Lists in SfdList_nom
  for (item in 1:n) {
    SListi   <- SfdList_nom[[item]]
    #  argument SfdList_nom already loaded with number of options M
    #  and nominal parameter matrix set up by mirt
    Mi      <- SListi$M
    # mirt parameter matrix for this item
    parmati <- SListi$Sfd
    #  compute exponentials of Mi linear basis functions
    expmati <- matrix(0,nfine,Mi)
    for (m in 1:Mi) expmati [,m] <- exp(parmati[m,1]*indfine + parmati[m,2])
    #  probability matrix
    probmati <- expmati/apply(expmati,1,sum)
    #  surprisal matrix
    surpmati <- -log(probmati)/log(Mi)
    #  approximate surprisal curves using spline basis
    Result <- smooth.basis(indfine, surpmati, Sbasis)
    Sfd <- Result$fd
    #  evaluate the first and second derivative matrices
    Dsurpmati  <- eval.fd(indfine, Sfd, 1)
    D2surpmati <- eval.fd(indfine, Sfd, 2)
    #  load these four matrices into SfdList_nom
    SListi$Pmatfine   <-   probmati
    SListi$Smatfine   <-   surpmati
    SListi$DSmatfine  <-  Dsurpmati
    SListi$D2Smatfine <- D2surpmati
    #  compute infoSurp mesh and length
    infovec           <- pracma::cumtrapz(indfine, sqrt(apply(Dsurpmati^2,1,sum)))
    SListi$infovec    <- infovec
    SListi$infoSurp   <- max(infovec)
    #  load the struct SList into List array
    SfdList_nom[[item]] <- SListi
  }
  return(SfdList_nom)
}
