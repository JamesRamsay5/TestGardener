\name{scoreDensity}
\alias{scoreDensity}
\title{
  Compute and plot a score density histogram and and curve.
}
\description{
  The tasks of function \code{theta.density()} and plotting the density are combined.  
  The score density is plotted both as a histogram and as a smooth curve.  
  All the score types may be plotted: sum scores, expected test scores, 
  percentile score index values, and locations on the test information or 
  scale curve.  The plot is output as a ggplot2 plot object, which is actually 
  plotted using the \code{print}
  command.
}
\usage{
  scoreDensity(scrvec, scrrng, ndensbasis=15, ttlstr=NULL, pltmax=0)
}
\arguments{
  \item{scrvec}{A vector of strictly increasing bin boundary values, with the 
  first at the lowest plotting value and the last at the upper boundary.  
  The number of bins in the histogram is one less than the number of bndry 
  values.}
  \item{scrrng}{A vector of length 2 containing lower and upper boundaries on 
  scores.}
  \item{ndensbasis}{The number of spline basis functions to be used to 
  represent the smooth density curve.}
  \item{ttlstr}{A string object used as a title for the plot.  
  Defaults to none.}
  \item{pltmax}{An upper limit on the vertical axis for plotting.  
  Defaults to the maximum curve value.}
}
\value{
A ggplot2 plot object \code{dens.plot} that can be displayed using command
\code{print(dens.plot)}.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{thetafun}},
  \code{\link{theta2arclen}},
  \code{\link{testscore}},
  \code{\link{theta.distn}}
}
\examples{
#  Example 1.  Display probability density curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList <- Quant_parList$WfdList
theta   <- Quant_parList$theta
#  plot the density for the score indices within interval c(0,100)
theta.scrrng <- c(0,100)
theta_int <- theta[0 < theta & theta < 100]
scoreDensity(theta_int, theta.scrrng)
#  plot the density for the sum score
sumscr <- Quant_dataList$scrvec
scoreDensity(sumscr, c(0,24), ttlstr="Sum score")
#  plot the density for the expected score
mu <- testscore(theta, WfdList, Quant_dataList$optList)
mu.scrrng <- c(0,24)
scoreDensity(mu, mu.scrrng, ttlstr="Expected sum score")
#  plot the density for the arclength score
ArclenResult <- theta2arclen(theta, WfdList, Quant_dataList$Wdim)
theta_al     <- ArclenResult$theta_al
arclen.scrrng <- c(0,ArclenResult$arclength)
scoreDensity(theta_al, arclen.scrrng, ttlstr="Arc length score")
#  Example 2.  Display probbility density curves for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}
