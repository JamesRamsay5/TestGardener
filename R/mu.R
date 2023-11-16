mu <- function(index, SfdList, scoreList){
  #  compute expected score mu
  
  # Last modified 9 November 2023 by Jim Ramsay
  
  N  <- length(index)
  n  <- length(SfdList)
  muvec <- rep(0,N)
  for (item in 1:n) {
    SListi <- SfdList[[item]]
    Sfdi   <- SListi$Sfd
    Mi     <- SListi$M
    if (Mi == 1){
      stop("Mi = 1.  Binary data should use Mi = 2.")
    } else {
      Smati <- eval.surp(index, Sfdi)
      Pmati <- exp(-Smati*log(Mi))
      scri  <- matrix(scoreList[[item]], N, Mi, byrow=TRUE)
      muvec <- muvec + apply(scri*Pmati,1,sum)
    }
  }
  return(muvec)
}
