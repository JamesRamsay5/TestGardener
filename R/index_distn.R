index_distn <- function(indexdens, logdensbasis,
                        pvec=c(0.05, 0.25, 0.50, 0.75, 0.95), nfine=101) {
  #  Compute cumulative density functions and probability density functions
  #  for the values in THETADENS using the basis object LOGDENSBASIS.

  # Last modified 12 December 2023 by Jim Ramsay
  
  #  check logdensfd
  
  if (!inherits(logdensbasis, "basisfd"))
    stop("Argument logdensbasis is not a functional basis object.")
  
  #  check indexdens
  
  if (!is.numeric(indexdens) && !inherits(indexdens,"matrix"))
    stop("Argument indexdens is neither a numeric vector or a matrix")
  if (inherits(indexdens,"matrix") && dim(indexdens)[2] > 1)
    stop("Argument is a matrix but has more than one column.")

  indrng  <- logdensbasis$rangeval # updated by Juan Li 2020-12-09
  indfine <- seq(indrng[1],indrng[2],len=nfine) # updated by Juan Li 2020-12-09

  #  set up basis and fd objects

  #  Estimate log density function and norming constant
  
  logdensfd    <- fda::fd(matrix(0,logdensbasis$nbasis,1), logdensbasis)
  rsList       <- TG_density.fd(indexdens, logdensfd, dbglev=0)
  logdensfd    <- rsList$Sfdobj
  C            <- rsList$C

  #  Compute probability and cumulative densities
  #  over a fine mesh of score values

  pdffine  <- exp(fda::eval.fd(indfine,logdensfd))/C
  cdffine  <- pracma::cumtrapz(indfine, pdffine)
  cdffine[cdffine < 0] <- 0
  cdffine[cdffine > 1] <- 1
  cdffine[1]     <- 0
  cdffine[nfine] <- 1
  
  #  convert to fd objects
  
  #  functional probability curves
  
  pdf_fd <- fda::smooth.basis(indfine, pdffine, logdensbasis)$fd

  #  compute cdf_fd using smooth.morph
  
  Snbasis <- 7
  Sbasis  <- fda::create.bspline.basis(indrng, Snbasis)
  Sfd     <- fda::fd(matrix(0,Snbasis,1),Sbasis)
  
  indexsort <- sort(indexdens[indexdens > indrng[1] & indexdens < indrng[2]])
  Nsort     <- length(indexsort)
  prbvec    <- (1:Nsort)/(Nsort+1)
  denscdf   <- unique(cdffine)
  indcdf    <- seq(indrng[1],indrng[2],len=length(denscdf))
  
  # return results
  
  return(list(pdf_fd=pdf_fd, cdffine=cdffine, pdffine=pdffine, 
              logdensfd=logdensfd, C=C, indcdf=indcdf, denscdf=denscdf))

}
