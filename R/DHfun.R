DHfun <- function(theta, WfdList, Umat) {

# Last modified 7 August 2023 by Jim Ramsay

  if (is.null(ncol(Umat)))
  {
    N <- 1
    n <- length(Umat)
  } else
  {
    if (ncol(Umat)==1)
    {
      N <- 1
      n <- length(Umat)
    } else
    {
      N <- nrow(Umat)
      n <- ncol(Umat)
    }
  }

  # loop through items to compute DH and D2H
  if (N == 1)
  {
    DH     <- 0
    D2H    <- 0
    Rveci  <- 0
    R2veci <- 0
  } else {
    DH     <- rep(0,N)
    D2H    <- rep(0,N)
    Rveci  <- rep(0,N)
    R2veci <- rep(0,N)
  }

  for (item in 1:n) {
    if (N == 1) {
      Uveci <- as.integer(Umat[item])
    } else {
      Uveci <- as.integer(Umat[,item])
    }

    if (!is.null(Uveci)) {
      #  extract the surprisal curves for this item
      WStri     <- WfdList[[item]]
      Wfdi      <- WStri$Wfd
      Mi        <- WStri$M
      #  evaluate surprisal curves at the score index values in theta
      DWmati    <- eval.surp(theta, Wfdi, 1)
      D2Wmati   <- eval.surp(theta, Wfdi, 2)
      #  Mi must be greater than 1, if not, abort
      if (Mi > 1) {
        #  select values of first and second derivatives of curve for the selected option
        if (N == 1) {
          Rveci  <-  DWmati[Uveci]
          R2veci <- D2Wmati[Uveci]
        } else {
          Wmati  <- rbind(DWmati,D2Wmati)
          for (j in 1:N)
          {
            Rveci[j]  <-  DWmati[j,Uveci[j]]
            R2veci[j] <- D2Wmati[j,Uveci[j]]
          }
        }
        # update fit derivative values
        DH  <- DH  +  Rveci
        D2H <- D2H + R2veci
      } else {
        stop("Mi = 1. Binary data should use Mi = 2.")
      }
    }
  }
  return(list(DH=DH, D2H=D2H))
}
