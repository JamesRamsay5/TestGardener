\name{theta.distn}
\alias{theta.distn}
\title{Compute score density}
\usage{
  theta.distn(thetadens, logdensbasis, nfine = 101)
}
\arguments{
  \item{thetadens}{A vector of score index, test score, or arc length values.  
  In the score index case, these are usually only the values in the interior of 
  the interval [0,100].}
 \item{logdensbasis}{A functional basis object for representing the log density 
 function.  The argument may also be a functional parameter object 
 (\code{fdPar}) or a functional basis object (\code{Wbasis}).}
  \item{nfine}{The number of values in a fine grid, default as 101.}
}
\description{
Computes the log density \code{logdensfd}, the norming constant \code{C} and the
values \code{denscdf} of the cumulative probability function over a set of score values.  
The score values may score index values, expected test score values, or locations on
the test information or scale curve.  The argument functional data object \code{logdensfd} 
should have a range that is appropriate for the score values being represented:
For score indices, [0,100], for expected test scores, the range of observed or 
expected scores; and for test information curve locations in the interval [0,\code{arclength}].
}
\value{
A named list containing:
  \describe{
    \item{logdensfd:}{A functional data object (\code{fd}) representing the log of the probability 
    function for input \code{theta}.} 
    \item{denscdf:}{A vector of the probability cumulative distribution function 
    values over the fine mesh.} 
    \item{C:}{The normalization constant for computing the probability 
    density function with the command \code{densityfd = exp(logdensfd)/C.}} 
    \item{densfine:}{The values of the probability density function over the fine mesh.}
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
WfdList   <- Quant_parList$WfdList
theta     <- Quant_parList$theta
N         <- length(theta)
#  Define the density for only interior theta values
inside    <- theta > 0 & theta < 100
logdensbasis    <- Quant_parList$logdensfd$basis
theta.distnList <- theta.distn(theta[inside], logdensbasis)
theta.densfine  <- theta.distnList$densfine
#  Add number of values at each boundary
thetazero <- length(theta[theta == 0])
thetacent <- length(theta[theta == 100])
#  Plot both the interior density and the boundary proportions
plot(indfine, theta.densfine, type="l", lwd=2, ylim=c(0,0.025), 
     xlab="Percentile Index", ylab="Density")
lines(c(  0,100), c(1/100,1/100), lty=2)
lines(c(  0,  0), c(0,thetazero/N), lty=1, lwd=4)
lines(c(100,100), c(0,thetacent/N), lty=1, lwd=4)
#  Second task:  Plot the test score density
optList <- Quant_dataList$optList
mu      <- testscore(theta, WfdList, optList)
#  Define the functional data object for score distribution over [0,24]
logdensmubasis <- fda::create.bspline.basis(c(0,24), 11)
logdensmufd    <- fd(matrix(0,11,1), logdensmubasis)
mu.distnList   <- theta.distn(mu, logdensmubasis)
mu.densfine    <- mu.distnList$densfine
scrfine        <- Quant_dataList$scrfine
#  Plot the test score density
plot(scrfine, mu.densfine, type="l", lwd=2, ylim=c(0,0.15), 
     xlab="Test Score", ylab="Expected Score Density")
}
