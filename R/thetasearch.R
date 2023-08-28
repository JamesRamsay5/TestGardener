thetasearch <- function(WfdList, U, theta, Hval, DHval, D2Hval, thetaind=1:N) {
  #  thetasearch examines the 101 points of each examinee's log likelihood 
  #  curve for (multiple minima.  if (the value in vector theta is more than 
  #  one away from the location of the minimum, theta is replaced by the
  #  index of the point and Hval, DHval, and D2Hval are revised.
  #  Integer changecount indicate the number of values changed.
  
  #  Last modified 14 August 2023
  
  # determine the number of rows in data matrix U and if (there is no subset
  #  specified in argument thetaind, search all of the rows of U.
  
  N <- nrow(U)
  n <- ncol(U)
  evalarg     <- seq(0,100,len=101)
  changeindex <- NULL
  for (j in thetaind) {
    #  identify the minima
    Hfine    <- Hcurve(WfdList, U[j,])
    Hfind    <- which.min(min(Hfine) == Hfine)
    Hminind  <- 0
    if (Hfind ==   1) Hminind <- Hfind
    if (Hfind == 101) Hminind <- Hfind
    if (Hminind == 0) {
      Hdiff    <- diff(Hfine)
      if (sum(Hdiff < 0) == 100){
        thetagrid <- 100
        Hvalgrid  <- Hfine[101]
      }
      if (sum(Hdiff > 0) == 100){
        thetagrid <- 0
        Hvalgrid  <- Hfine[1]
      }
      if (sum(Hdiff < 0) < 100 && sum(Hdiff > 0) < 100) {
        #  normal case:  find the minimum with the lowest value
        thetagrid <- evalarg[Hfine == min(Hfine)]
        Hvalgrid  <- Hfine[thetagrid+1]
      }
      #  now check that the input value of theta is close to thetagrid
      change <- theta[j] - thetagrid
      if (abs(change) >= 0.5) {
        #  the input theta and the grid value are different.
        #  replace the input value with the grid value
        theta[j]  <- thetagrid
        Hval[j]   <- Hvalgrid
        DHval[j]  <- 0
        D2Hval[j] <- 0
        for (item in 1:n) {
          WListi      <- WfdList[[item]]
          DWmatfinei  <- WListi$DWmatfine
          DHval[j]    <-  DHval[j] +  DWmatfinei[thetagrid+1,U[j,item]]
          D2Wmatfinei <- WListi$D2Wmatfine
          D2Hval[j]   <- D2Hval[j] + D2Wmatfinei[thetagrid+1,U[j,item]]
        }
        changeindex <- c(changeindex,j)
      }
    }
  }
  
  return(list(theta=theta, Hval=Hval, DHval=DHval, D2Hval=D2Hval, 
              changeindex=changeindex))
  
}


