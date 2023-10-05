\name{TestInfo_svd}
\alias{TestInfo_svd}
\title{
  Image of the Test Tnformation Curve in 2 or 3 Dimensions
}
\description{
  The test information curve is the trajectory of joint variation of all 
  the surprisal curves within the ambient space of dimension the total number 
  of curves.  But usually a very high percent of the shape variation in the 
  curve can be represented in either two or three dimensions using the singular
  value decomposition of a matrix of total curve values over a fine mesh.  
  The resulting approximation is converted to a set of surprisal curve values.
}
\usage{
  TestInfo_svd(scrfine, WfdList, itemindex=1:n, nharm=2)
}
\arguments{
  \item{scrfine}{A fine mesh of values over which the image is plotted.  This 
  is usually either the score index theta or the test arc length.}
  \item{WfdList}{A list vector of length n, the number of test items.  
  Each list in the vector contains values of the surprisal curves for that
  item.}
  \item{itemindex}{A vector of item indices to be used in the approximation.}
  \item{nharm}{The number of dimension in the approximation, usually either 
  two or three.}
}
\value{
  The approximation is returned as a surprisal functional data object, and so
  are the percentages of the total variation fit by each dimension in the
  approximation.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
