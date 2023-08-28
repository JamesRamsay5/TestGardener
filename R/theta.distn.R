theta.distn <- function(thetadens, logdensbasis,
                        pvec=c(0.05, 0.25, 0.50, 0.75, 0.95), nfine=101) {
  #  Compute cumulative density functions and probability density functions
  #  for the values in THETADENS using the basis object LOGDENSBASIS.

  # Last modified 7 March 2023 by Jim Ramsay
  
  #  check logdensfd
  
  if (!inherits(logdensbasis, "basisfd"))
    stop("Argument logdensbasis is not a functional basis object.")
  
  #  check thetadens
  
  if (!is.numeric(thetadens) && !inherits(thetadens,"matrix"))
    stop("Argument thetadens is neither a numeric vector or a matrix")
  if (inherits(thetadens,"matrix") && dim(thetadens)[2] > 1)
    stop("Argument is a matrix but has more than one column.")

  indrng  <- logdensbasis$rangeval # updated by Juan Li 2020-12-09
  indfine <- seq(indrng[1],indrng[2],len=nfine) # updated by Juan Li 2020-12-09

  #  set up basis and fdPar objects

  #  Estimate log density function and norming constant
  
  logdensfd    <- fda::fd(matrix(0,logdensbasis$nbasis,1), logdensbasis)
  logdensfdPar <- fda::fdPar(logdensfd)
  rsList       <- TG_density.fd(thetadens, logdensfdPar, dbglev=0)
  logdensfd    <- rsList$Wfdobj
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
  
  pdf_fd <- smooth.basis(indfine, pdffine, logdensbasis)$fd
  # cdf_fd <- smooth.basis(indfine, cdffine, logdensbasis)$fd

  #  compute cdf_fd using smooth.morph
  
  Wnbasis <- 7
  Wbasis  <- create.bspline.basis(indrng, Wnbasis)
  Wfd     <- fd(matrix(0,Wnbasis,1),Wbasis)
  WfdPar  <- fdPar(Wfd)
  
  thetasort <- sort(thetadens[thetadens > indrng[1] & thetadens < indrng[2]])
  Nsort     <- length(thetasort)
  prbvec    <- (1:Nsort)/(Nsort+1)
  
  # result    <- smooth.morph(thetasort, prbvec, c(0,1), WfdPar)
  # ysmth     <- result$ysmth 
  # thetafine <- c(indrng[1], thetasort, indrng[2])
  # yfine     <- c(0, ysmth, 1   )
  
  denscdf   <- unique(cdffine)
  indcdf    <- seq(indrng[1],indrng[2],len=length(denscdf))
  
  # return results
  
  return(list(pdf_fd=pdf_fd, cdffine=cdffine, pdffine=pdffine, 
              logdensfd=logdensfd, C=C, indcdf=indcdf, denscdf=denscdf))

}
