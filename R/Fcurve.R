Fcurve <- function(SfdList, chcevec, nderiv=0) {
  #  FCURVE computes 101 values one of the first three derivatives of the  
  #. negative log likelihood for a set of examinees over [the width of 
  #. the interval [0,100] at 101 equally spaced points..
  #. Arguments:
  #  SfdList   ... A list vector containing descriptions of each item
  #. The curve is defined for each rater by the choices made.
  #  chcevec         ... The N by n matrix of choice indices
  #. nderiv    ... The order of derivative selected from 0, 1 or 2.
  
  #  Last updated 3 November 2023
  
  #. Deal with a single choice vector that is numeric
  if (!is.matrix(chcevec)) chcevec <- t(as.matrix(chcevec))
  nitem <- ncol(chcevec)
  
  #  loop through items to compute negative log likelihood values in F
  
  Ffine <- matrix(0,101,1)
  for (item in 1:nitem) {
    chceveci <- chcevec[,item]
    if (!is.null(chceveci)) {
      SListi    <- SfdList[[item]]
      if (nderiv == 0) {
        Smati <- SListi$Smatfine
        Ffine <- Ffine + Smati[,chceveci]
      }
      if (nderiv == 1) {
        DSmati <- SListi$DSmatfine
        Ffine  <- Ffine + DSmati[,chceveci]
      }
      if (nderiv == 2) {
        D2Smati <- SListi$D2Smatfine
        Ffine   <- Ffine + D2Smati[,chceveci]
      }
    }
  }
  
  return(Ffine)
  
}
