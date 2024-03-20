DFfun <- function(index, SfdList, chcemat) {

# Last modified 19 December 2023 by Jim Ramsay

  if (is.null(ncol(chcemat)))
  {
    N <- 1
    n <- length(chcemat)
  } else {
    if (ncol(chcemat)==1)
    {
      N <- 1
      n <- length(chcemat)
    } else {
      N <- nrow(chcemat)
      n <- ncol(chcemat)
    }
  }

  # loop through items to compute DF and D2F
  if (N == 1) {
    DF     <- 0
    D2F    <- 0
    Rveci  <- 0
    R2veci <- 0
  } else {
    DF     <- rep(0,N)
    D2F    <- rep(0,N)
    Rveci  <- rep(0,N)
    R2veci <- rep(0,N)
  }

  for (item in 1:n) {
    if (N == 1) {
      chceveci <- as.integer(chcemat[item])
    } else {
      chceveci <- as.integer(chcemat[,item])
    }

    # if (!is.null(chceveci)) {
      #  extract the surprisal curves for this item
      SStri     <- SfdList[[item]]
      Sfdi      <- SStri$Sfd
      Mi        <- SStri$M
      Zmati     <- SStri$Zmat
      #  evaluate surprisal curves at the score index values in index
      DSmati    <- eval.surp(index, Sfdi, Zmati, 1)
      D2Smati   <- eval.surp(index, Sfdi, Zmati, 2)
      #  Mi must be greater than 1, if not, abort
      # if (Mi > 1) {
        #  select values of first and second derivatives of curve for the selected option
        if (N == 1) {
          Rveci  <-  DSmati[chceveci]
          R2veci <- D2Smati[chceveci]
        } else {
          Smati  <- rbind(DSmati,D2Smati)
          for (j in 1:N) {
            Rveci[j]  <-  DSmati[j,chceveci[j]]
            R2veci[j] <- D2Smati[j,chceveci[j]]
          }
        }
        # update fit derivative values
        DF  <- DF  +  Rveci
        D2F <- D2F + R2veci
      # } else {
      #   stop("Mi = 1. Binary data should use Mi = 2.")
      # }
    # }
  }
  return(list(DF=DF, D2F=D2F))
}
