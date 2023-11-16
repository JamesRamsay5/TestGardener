entropies <- function(index, m, n, chcemat, noption) {
#  Compute two single entropies, joint entropy, and mutual entropy
#  arguments
#  index   ... vector of score index values
#  m       ... index of first item
#  n       ... index of second item
#  chcemat       ... data matrix containing choice indices
#  noption ... number of options in an item
  
#  Last modified 2 November 2023 by Jim
  
  nobs  <- length(index)
  nobsinv <- 1/nobs
  #  single entropies
  #  single entropy m
  Hmprob <- matrix(0,noption[m],1)
  for (j in 1:nobs) {
    mi <- chcemat[j,m]
    Hmprob[mi] <- Hmprob[mi] + nobsinv
  }
  Hm <- 0
  for (mi in 1:noption[m]) {
    Pm <- Hmprob[mi]
    Hm <- Hm - Pm*log(Pm)
  }
  #  single entrop n
  Hnprob <- matrix(0,noption[n],1)
  for (j in 1:nobs) {
    ni <- chcemat[j,n]
    Hnprob[ni] <- Hnprob[ni] + nobsinv
  }
  Hn <- 0
  for (ni in 1:noption[n]) {
    Pn <- Hnprob[ni]
    Hn <- Hn - Pn*log(Pn)
  }
  #  joint entropy mn
  Hjointprob <- matrix(0,noption[m],noption[n])
  for (j in 1:nobs) {
    mi <- chcemat[j,m]
    ni <- chcemat[j,n]
    Hjointprob[mi,ni] <- Hjointprob[mi,ni] + nobsinv
  }
  Hjoint <- 0
  for (mi in 1:noption[m]) {
    for (ni in 1:noption[n]) {
      Pmn <- Hjointprob[mi,ni]
      if (Pmn > 0) {
        Hjoint <- Hjoint - Pmn*log(Pmn)
      }
    }
  }
  #  mutual entropy 
  Hmutual <- Hm + Hn - Hjoint
  
  return(list(Hm=Hm, Hn=Hn, Hjoint=Hjoint, Hmutual=Hmutual))
  
}