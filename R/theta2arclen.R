theta2arclen <- function(theta, WfdList, Wdim) {
  
  # Last modified 8 February 2021 by Jim Ramsay

  #   set up a fine mesh over score index range, which is assumed to be [0,100]
  nfine   <- 1001
  indfine <- seq(0,100,length.out=nfine)
  
  #  compute W derivatives over fine mesh
  DWfine  <- matrix(0,nfine,Wdim)
  n       <- length(WfdList)
  m2      <- 0
  for (i in 1:n)
  {
    WListi <- WfdList[[i]]
    Wfdi  <- WListi$Wfd
    Mi    <- WListi$M
    m1    <- m2 + 1
    m2    <- m2 + Mi
    DWfinei <- eval.surp(indfine, Wfdi, 1)
    DWfine[,m1:m2] <- DWfinei
  }
  # accumulate arclength values
  delta         <- 0.1
  arclengthfine <- delta * as.numeric(pracma::cumtrapz(sqrt(rowSums(DWfine^2))))
  # add zero to complete the mesh and arc length values
  arclength     <- max(arclengthfine)
  theta_al      <- pracma::interp1(indfine, arclengthfine, theta)
  Qvec_al <- pracma::interp1(indfine, arclengthfine, c(5,25,50,75,95))
  arclengthfine <- arclengthfine[seq(1,1001,10)]
  return(list(theta_al=theta_al, arclength=arclength, arclengthfine=arclengthfine, 
              Qvec_al=Qvec_al))
}