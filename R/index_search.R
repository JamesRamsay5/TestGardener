index_search <- function(SfdList, chcemat, index, Fval, DFval, D2Fval, 
                         indexind=1:N) {
  #  indexs_earch examines the 101 points of each examinee's log likelihood 
  #  curve for (multiple minima.  if (the value in vector index is more than 
  #  one away from the location of the minimum, index is replaced by the
  #  index of the point and Fval, DFval, and D2Fval are revised.
  #  Integer changecount indicate the number of values changed.
  
  #  Last modified 21 November 2023
  
  # determine the number of rows in data matrix chcemat and if (there is no 
  #  subset specified in argument indexind, search all of the rows of chcemat.
  
  N <- nrow(chcemat)
  n <- ncol(chcemat)
  evalarg     <- seq(0,100,len=101)
  changeindex <- NULL
  for (j in indexind) {
    # print(j)
    #  identify the minima
    chcevecj <- chcemat[j,]
    Ffine    <- Fcurve(SfdList, chcevecj)
    Ffind    <- which.min(min(Ffine) == Ffine)
    Fminind  <- 0
    if (Ffind ==   1) Fminind <- Ffind
    if (Ffind == 101) Fminind <- Ffind
    if (Fminind == 0) {
      Fdiff    <- diff(Ffine)
      if (sum(Fdiff < 0) == 100){
        indexgrid <- 100
        Fvalgrid  <- Ffine[101]
      }
      if (sum(Fdiff > 0) == 100){
        indexgrid <- 0
        Fvalgrid  <- Ffine[1]
      }
      if (sum(Fdiff < 0) < 100 && sum(Fdiff > 0) < 100) {
        #  normal case:  find the minimum with the lowest value
        indexgrid <- evalarg[Ffine == min(Ffine)]
        Fvalgrid  <- Ffine[indexgrid+1]
      }
      #  now check that the input value of index is close to indexgrid
      change <- index[j] - indexgrid
      if (abs(change) >= 0.5) {
        #  the input index and the grid value are different.
        #  replace the input value with the grid value
        index[j]  <- indexgrid
        Fval[j]   <- Fvalgrid
        DFval[j]  <- 0
        D2Fval[j] <- 0
        for (item in 1:n) {
          SListi      <- SfdList[[item]]
          DSmatfinei  <- SListi$DSmatfine
          DFval[j]    <- DFval[j] +  DSmatfinei[indexgrid+1,chcemat[j,item]]
          D2Smatfinei <- SListi$D2Smatfine
          D2Fval[j]   <- D2Fval[j] + D2Smatfinei[indexgrid+1,chcemat[j,item]]
        }
        changeindex <- c(changeindex,j)
      }
    }
  }
  
  return(list(index=index, Fval=Fval, DFval=DFval, D2Fval=D2Fval, 
              changeindex=changeindex))
  
}


