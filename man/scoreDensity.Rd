\name{scoreDensity}
\alias{scoreDensity}
\title{
  Compute and plot a score density histogram and and curve.
}
\description{
  The tasks of function \code{index.density()} and plotting the density are combined.  
  The score density is plotted both as a histogram and as a smooth curve.  
  All the score types may be plotted: sum scores, expected test scores, 
  percentile score index values, and locations on the test information or 
  scale curve.  The plot is output as a ggplot2 plot object, which is actually 
  plotted using the \code{print}
  command.
}
\usage{
  scoreDensity(scrvec, scrrng=c(0,100), ndensbasis=15, ttlstr=NULL, pltmax=0)
}
\arguments{
  \item{scrvec}{A vector of strictly increasing bin boundary values, with the 
    first at the lowest plotting value and the last at the upper boundary.  
    The number of bins in the histogram is one less than the number of bndry 
    values.}
  \item{scrrng}{A vector of length 2 containing lower and upper boundaries on 
    scores, which defaults to c(0,100).}
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
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{index_fun}},
  \code{\link{index2info}},
  \code{\link{mu}},
  \code{\link{index_distn}}
}
\examples{
#  Example 1.  Display probability density curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
SfdList <- Quant_13B_problem_parmList$SfdList
index   <- Quant_13B_problem_parmList$index
Qvec    <- Quant_13B_problem_parmList$Qvec
#  plot the density for the score indices within interval c(0,100)
index_int <- index[0 < index & index < 100]
oldpar <- par(no.readonly=TRUE)
scoreDensity(index_int)
par(oldpar)
}
