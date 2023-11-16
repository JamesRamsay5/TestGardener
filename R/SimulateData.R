SimulateData <- function(nsim, indfine, denscdf, SfdList) {

  # Last modified 31 October 2023 by Jim Ramsay
  
  n       <- length(SfdList)
  noption <- rep(0,n)
  for (i in 1:n) noption[i] <- SfdList[[i]]$M
  ndens   <- length(denscdf)
  Usim    <- matrix(0,nsim,n)
  nj2 <- 0
  for (j in 2:ndens) {
    nj1 <- nj2 + 1
    if (j < ndens) {
      nthetaj <- round((denscdf[j] - denscdf[j-1])*nsim,0)
      nj2 <- nj2 + nthetaj 
    }  else {
      nj2 <- nsim
      nthetaj <- nj2-nj1+1
    }
    thetaj <- (indfine[j] + indfine[j-1])/2
    # print(c(j,thetaj,nthetaj,nj2))
    for (i in 1:n) {
      WListi <- SfdList[[i]]
      Mi <- WListi$M
      ind <- 1:Mi
      surpij <- eval.surp(thetaj, WListi$Sfd)
      probij <- Mi^(-surpij)
      for (k in 1:nthetaj) {
        choiceij <- rmultinom(1, Mi, probij)
        Usim[nj1+k-1,i] <- sum(ind*choiceij)
      }
    }
  }
  return(Usim)
}
