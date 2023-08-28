testscore <- function(theta, WfdList, optList){
  #  compute expected score mu
  
  # Last modified 7 August 2023 by Jim Ramsay
  
  N  <- length(theta)
  n  <- length(WfdList)
  mu <- rep(0,N)
  for (item in 1:n) {
    SListi <- WfdList[[item]]
    Sfdi   <- SListi$Wfd
    Mi     <- SListi$M
    if (Mi == 1){
      stop("Mi = 1.  Binary data should use Mi = 2.")
    } else {
      Smati <- eval.surp(theta, Sfdi)
      Pmati <- exp(-Smati*log(Mi))
      scri  <- matrix(optList$optScr[[item]], N, Mi, byrow=TRUE)
      mu    <- mu + apply(scri*Pmati,1,sum)
    }
  }
  return(mu)
}
