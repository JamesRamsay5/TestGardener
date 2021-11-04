theta.distn <- function(thetadens, logdensbasis, nfine = 101) {

# Last modified 8 February 2021 by Jinm Ramnsay

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

  densfine <- exp(fda::eval.fd(indfine,logdensfd))/C
  delta    <- indfine[2] - indfine[1]
  denscdf  <- delta*pracma::cumtrapz(densfine)
  denscdf[denscdf < 0] <- 0
  denscdf[denscdf > 1] <- 1
  denscdf[1]     <- 0
  denscdf[nfine] <- 1

  # return results
  
  return(list(logdensfd=logdensfd, denscdf=denscdf, C=C, densfine=densfine))

}
