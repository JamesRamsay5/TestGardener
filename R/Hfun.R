
Hfun <- function(theta, WfdList, Umat) {
	
# Last modified 4 December 2020 by Jim Ramsay

  if (is.null(ncol(Umat))) {
    N <- 1
    n <- length(Umat)
  } else {
    if (ncol(Umat)==1) {
      N <- 1
      n <- length(Umat)
    } else {
      N <- nrow(Umat)
      n <- ncol(Umat)
    }
  }
  
  #  set up vectors to contain fitand curve values
  
  if (N == 1) {
    H     <- 0
    Wveci <- 0
  } else {
    H     <- rep(0,N)
    Wveci <- rep(0,N)
  }

  # loop through items to compute negative log likelihood values in H
  
  for (item in 1:n) {
    #  be sure that U contains only integers
    if (N == 1) {
      Uveci <- as.integer(Umat[item])
    } else {
      Uveci <- as.integer(Umat[,item])
    }
    #  Now compute increment to fit values for this item
    #  provided Uveci is not NULL
    if (!is.null(Uveci)) {
      #  extract the surprisal curves for this item
      WStri     <- WfdList[[item]]
      Wfdi      <- WStri$Wfd
      Mi        <- WStri$M
      #  evaluate surprisal curves at the score index values in theta
      Wmati     <- eval.surp(theta, Wfdi)
      #  Mi must be greater than 1, if not, abort
      if (Mi > 1) {
        #  select values of curve for the selected option
        if (N == 1) {
          Wveci <- Wmati[Uveci]
        } else {
          for (j in 1:N)
          {
            if (Uveci[j] > Mi)
            {
              stop(paste("Item: ", item, " Uveci(",j,") > Mi",sep = ""))
            }
            Wveci[j] <- Wmati[j,Uveci[j]]
          }
        }
        # update fit values
        H <- H + Wveci
      } else {
        stop("Mi not greater than 1. Binary data should use Mi = 2.")
      }
    }
  }
  return(H)
}
