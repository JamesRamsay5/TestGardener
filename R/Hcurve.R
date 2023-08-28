Hcurve <- function(WfdList, U, nderiv=0) {
  #  HCURVE computes 101 values one of the first three derivatives of the  
  #. negative log likelihood for a set of examinees over [the width of 
  #. the interval [0,100] at 101 equally spaced points..
  #. Arguments:
  #  WfdList   ... A list vector containing descriptions of each item
  #. The curve is defined for each rater by the choices made.
  #  U         ... The N by n matrix of choice indices
  #. nderiv    ... The order of derivative selected from 0, 1 or 2.
  
  #  Last updated 15 August 2023
  
  #. Deal with a single choice vector that is numeric
  if (!is.matrix(U)) U <- t(as.matrix(U))
  nitem   <- ncol(U)
  
  #  loop through items to compute negative log likelihood values in H
  
  Nperson <- nrow(U)
  Hfine <- matrix(0,101,Nperson)
  for (item in 1:nitem) {
    Uveci <- U[,item]
    if (!is.null(Uveci)) {
      WListi    <- WfdList[[item]]
      if (nderiv == 0) {
        Wmati <- WListi$Wmatfine
        Hfine <- Hfine + Wmati[,Uveci]
      }
      if (nderiv == 1) {
        DWmati <- WListi$DWmatfine
        Hfine  <- Hfine + DWmati[,Uveci]
      }
      if (nderiv == 2) {
        D2Wmati <- WListi$D2Wmatfine
        Hfine   <- Hfine + D2Wmati[,Uveci]
      }
    }
  }
  
  return(Hfine)
  
}
