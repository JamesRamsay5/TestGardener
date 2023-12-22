Ffun <- function(index, SfdList, chcemat) {
	
# Last modified 19 December 2023 by Jim Ramsay

  if (is.null(ncol(chcemat))) {
    N <- 1
    n <- length(chcemat)
  } else {
    if (ncol(chcemat)==1) {
      N <- 1
      n <- length(chcemat)
    } else {
      N <- nrow(chcemat)
      n <- ncol(chcemat)
    }
  }
  
  #  set up vectors to contain fitand curve values
  
  if (N == 1) {
    F     <- 0
    Sveci <- 0
  } else {
    F     <- rep(0,N)
    Sveci <- rep(0,N)
  }

  # loop through items to compute negative log likelihood values in F
  
  for (item in 1:n) {
    #  be sure that chcemat contains only integers
    if (N == 1) {
      chceveci <- as.integer(chcemat[item])
    } else {
      chceveci <- as.integer(chcemat[,item])
    }
    #  Now compute increment to fit values for this item
    #  provided chceveci is not NULL
    if (!is.null(chceveci)) {
      #  extract the surprisal curves for this item
      SStri     <- SfdList[[item]]
      Sfdi      <- SStri$Sfd
      Mi        <- SStri$M
      Zmati     <- SStri$Zmat
      #  evaluate surprisal curves at the score index values in theta
      Smati     <- eval.surp(index, Sfdi, Zmati)
      #  Mi must be greater than 1, if not, abort
      if (Mi > 1) {
        #  select values of curve for the selected option
        if (N == 1) {
          Sveci <- Smati[chceveci]
        } else {
          for (j in 1:N) {
            # if (chceveci[j] > Mi) {
            #   stop(paste("Item: ", item, " chceveci(",j,") > Mi",sep = ""))
            # }
            Sveci[j] <- Smati[j,chceveci[j]]
          }
        }
        # update fit values
        F <-  + Sveci
      } else {
        stop("Mi not greater than 1. Binary data should use Mi = 2.")
      }
    }
  }
  return(F)
}
