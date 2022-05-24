density_plot <- function(scrvec, scrrng, Qvec, xlabstr=NULL, titlestr=NULL, 
                         scrnbasis=15, nfine=101) {
  #  DENSITY_PLOT plots the probability density function of a set of 
  #  score values that are not at the score boundaries as a smooth
  #  curve, and also plots the proportions of score values at both
  #  boundaries as points.  
  #  The score values are typically either the values of the score index 
  #  values theta or the arclength or information score values.  
  #  Arguments:
  #  SCRVEC   ... A vector of N score values
  #  SCRNG    ... A vector of length 2 containing boundary values
  #  QVEC     ... A vector of length 5 containing the score values 
  #               corresponding to the marker percentages 
  #               5, 25, 50, 75 and 95.
  #  XLABSTR  ... Label for abscissa
  #  TITLESTR ... Label for plot
  #  SCRNASIS ... The number of spline basis functions used for representing
  #               the smooth density function.
  #  NFINE    ... The number of plotting points.
  
  #  Last modified 18 May 2022 by Jim Ramsay
  
  #  set default values
  #  get the score values not on the boundaries
  N <- length(scrvec)
  scrdens <- scrvec[scrrng[1] < scrvec & scrvec < scrrng[2]]
  #  set up the basis object for the spline function representing
  #  the density function
  logdensbasis <- create.bspline.basis(scrrng, scrnbasis)    
  #  compute the values of the density function at a 
  densResults  <- theta.distn(scrdens, logdensbasis)
  logdensfd    <- densResults$logdensfd
  C            <- densResults$C
  densfine     <- densResults$pdffine
  #  plot the density
  scrfine      <- seq(scrrng[1],scrrng[2],len=nfine)
  N_max <- sum(scrvec == scrrng[2])
  N_min <- sum(scrvec == scrrng[1])
  pmax  <- max(c(N_min/N, max(densfine),N_max/N)) 
  plot(scrfine, densfine, type="l", lwd=2, ylim=c(0,pmax),
       xlab=xlabstr, ylab="Density", main=titlestr)
  for (k in 1:5) lines(c(Qvec[k],Qvec[k]), c(0,pmax), lty=2)
  points(scrrng[1], N_min/N, pch="o", lwd=2)
  points(scrrng[2], N_max/N, pch="o", lwd=2)
}