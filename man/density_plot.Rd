\name{density_plot}
\alias{density_plot}
\title{
  Plot the probability density function for a set of test scores 
}
\description{
Plots the probability density function of a set of 
score values that are not at the score boundaries as a smooth
curve, and also plots the proportions of score values at both
boundaries as points.  
The score values are typically either the values of the score index 
values theta or the arclength or information score values.  
}
\usage{
density_plot(scrvec, scrrng, Qvec, xlabstr=NULL, titlestr=NULL, 
                         scrnbasis=15, nfine=101)
}
\arguments{
  \item{scrvec}{A vector of N score values}
  \item{scrrng}{A vector of length 2 containing boundary values}
  \item{Qvec}{A vector of length 5 containing the score values 
  corresponding to the marker percentages 5, 25, 50, 75 and 95.}
  \item{xlabstr}{Label for abscissa}
  \item{titlestr}{Label for plot}
  \item{scrnbasis}{The number of spline basis functions used 
  for representing the smooth density function}
  \item{nfine}{Number of plotting points}
}
\value{
  A plot object.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{
Juan Li and James Ramsay
}
\seealso{
  \code{\link{scoreDensity}}
}
\examples{
#  Example 1.  Display probability density curve for the 
#  score index values for the short SweSAT multiple choice test with 
#  24 items and 1000 examinees
theta   <- Quantshort_parList$theta
Qvec    <- Quantshort_parList$Qvec
#  plot the density for the score indices within interval c(0,100)
density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
             titlestr="SweSAT 13B Theta Density",  
             scrnbasis=11, nfine=101)
#  Example 2.  Display probability density curve for the 
#  arc length or information values
theta_al  <- Quantshort_infoList$theta_al
Qvec_al   <- Quantshort_infoList$Qvec_al
arclength <- Quantshort_infoList$arclength
#  plot the density for the score indices within interval c(0,arclength)
density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Score index", 
             titlestr="SweSAT 13B Theta Density",  
             scrnbasis=11, nfine=101)
}
