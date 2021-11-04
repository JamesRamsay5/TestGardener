\name{ArcLength.plot}
\alias{ArcLength.plot}
\title{
Plot arc length as a function of score index \code{theta}.
}
\description{
Arc length is the distance along the space curved traced out as score index 
\code{theta} increases from 0 to 100.  It is measured in bits and is remains 
unchanged if the score index continuum is modified.}
\usage{
  ArcLength.plot(arclength, alfine, titlestr=NULL)
}
\arguments{
  \item{arclength}{This is a positive real number indicating the total length of the 
  space curve.  It is expressed in terms of numbers of bits.}
  \item{alfine}{A vector of length 101 containing equally-spaced arc-length distances  
  along the test information curve.}
  \item{titlestr}{A string for the title of the data.}
}
\value{
A \code{gg} or \code{ggplot} object defining the plot of arclength along the test information 
curve as a function of the score index \code{theta}.  This is displayed by the \code{print} command. 
The plot is automatically displayed as a side value even if no return object is specified in 
the calling statement.
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
  \code{\link{theta2arclen}}
}
\examples{
#  Example 1.  Display the arc length curve for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
indfine       <- seq(0,100,len=101)
WfdList       <- Quant_parList$WfdList
theta         <- Quant_parList$theta
alfine        <- Quant_parList$alfine
arclenList    <- theta2arclen(theta, WfdList, Quant_dataList$Wdim)
names(arclenList)
ArcLength.plot(arclenList$arclength, alfine)
#  Example 2.  Display the arc length curve for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}

