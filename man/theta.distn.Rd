\name{theta.distn}
\alias{theta.distn}
\title{Compute score density}
\usage{
  theta.distn(thetadens, logdensbasis,  
              pvec=c(0.05, 0.25, 0.50, 0.75, 0.95), nfine = 101)
}
\arguments{
  \item{thetadens}{A vector of score index, test score, or arc length values.  
    In the score index case, these are usually only the values in the interior of 
    the interval [0,100].}
  \item{logdensbasis}{A functional basis object for representing the log density 
    function.  The argument may also be a functional parameter object 
    (\code{fdPar}) or a functional basis object (\code{Wbasis}).}
  \item{pvec}{A vector length NL containing the marker percentages.}
  \item{nfine}{The number of values in a fine grid, default as 101.}
}
\description{
Computes the cumulated density for distribution function,
the probability density function, and the log probability density function
as fd objects by spline smoothing of the score values \code{thetadens} using 
the basis object \code{logdensbasis}.  The norming constant \code{C} is also 
output.

The score values may score index values \code{theta}, expected test score 
values \code{mu}, or arc length locations on the test information or scale 
curve.  The argument functional data object \code{logdensfd} should have a 
range that is appropriate for the score values being represented:
For score indices, [0,100], for expected test scores, the range of observed or 
expected scores; and for test information curve locations in the interval 
[0,\code{arclength}].
}
\value{
A named list containing:
  \describe{
    \item{pdf_fd:}{An fd object for the probability density function
      values over the fine mesh.} 
    \item{cdffine:}{A vector of cumulative probability values beginning
    with zero and ending with 1. It must not have ties.}
    \item{pdffine:}{A vector of probability values.}
    \item{logdensfd:}{A functional data object (\code{fd}) representing 
      the log of the probability function for input \code{theta}.} 
    \item{C:}{The normalization constant for computing the probability 
      density function with the command \code{densityfd = exp(logdensfd)/C.}} 
    \item{denscdf:}{A set of unique values of the cumulative probability
      function defined over an equally spaced mesh of score index values
      of the same length as \code{denscdf}.}
    \item{indcdf:}{A vector of values within [0,100] corresponding to the 
    values in \code{denscdf}.}
  }
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{thetafun}},
  \code{\link{theta2arclen}},
  \code{\link{testscore}},
  \code{\link{scoreDensity}}
}
\examples{
#  Example 1.  Display the item power curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  Assemble information for estimating theta density
indfine   <- seq(0,100,len=101)
WfdList   <- Quant_13B_problem_parList$WfdList
theta     <- Quant_13B_problem_parList$theta
N         <- length(theta)
#  Define the density for only interior theta values
inside    <- theta > 0 & theta < 100
logdensbasis    <- Quant_13B_problem_parList$logdensfd$basis
theta.distnList <- theta.distn(theta[inside], logdensbasis)
theta.densfine  <- eval.fd(indfine, theta.distnList$pdf_fd)
#  Add number of values at each boundary
thetazero       <- length(theta[theta == 0])
thetacent       <- length(theta[theta == 100])
#  Plot both the interior density and the boundary proportions
oldpar <- par(no.readonly=TRUE)
on.exit(oldpar)
plot(indfine, theta.densfine, type="l", lwd=2, ylim=c(0,0.025), 
     xlab="Percentile Index", ylab="Density")
lines(c(  0,100), c(1/100,1/100), lty=2)
lines(c(  0,  0), c(0,thetazero/N), lty=1, lwd=4)
lines(c(100,100), c(0,thetacent/N), lty=1, lwd=4)
}
