\name{theta2arclen}
\alias{theta2arclen}
\title{
  Compute results using arc length or information as the abscissa.
}
\description{The one-dimensional psychometric model defines a space curve within the vector space
  defined by the total collection of option surprisal curves.  This curve is a valuable
  resource since positions along the curve are defined in bits and positions on the
  curve are subject to the same strict properties that apply to physical measurements.

  Function \code{theta2arclen} is required to convert
  objects defined over the score index continuum \code{c(0,100)} to
  the same objects over the arc length continuum \code{c(0,arclength)}, and also
  vice versa.  Since the arc length or information continuum is along
  a space curve that is invariant under strictly monotone transformations
  of the score index \code{theta}, and is also a metric, it is an ideal 
  choice for the abscissa in all plots.}
\usage{
   theta2arclen(theta, Qvec, WfdList, binctr, itemindex=1:n, plotrng=c(0,100), shortwrd)
}
\arguments{
  \item{theta}{A vector of score index, test score, or arc length values, one for
    each examinee or respondent.}
  \item{Qvec}{A vector of locations of the five marker percentages.}
  \item{WfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  Its length
    is equal to the number of items in the test or questions in the scale.  
    Each member of \code{WfdList} is a named list containing information computed during the analysis.
    These named lists contain these objects:
    \describe{
      \item{Wfd:}{A functional data object containing the \code{M} surprisal curves
        for a question.}
      \item{M:}{The number of options.}
      \item{Pbin:}{A matrix containing proportions at each bin.}
      \item{Wbin:}{A matrix containing surprisal values at each bin.}
      \item{Pmatfine:}{A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:}{A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:}{A matrix of the values of the first derivative of 
        surprisal curves over fine mesh.}
      \item{D2Wmatfine:}{A matrix of the values of the second derivative of 
        surprisal curves over fine mesh.}
    }
  }
  \item{binctr}{A vector of locations of the bin centers.}
  \item{itemindex}{A vector containing the indices of the items to be used.}
  \item{plotrng}{A vector of length 2 containing the starting score index
    and end score index values of the range to be plotted.}
  \item{shortwrd}{If TRUE only vectors arclength and arclengtvec are returned in order
         to speed up the computation within cycles in function Analyze where only these 
         objects are required.  The default is FALSE.}
}
\value{A named list object containing these results of the analysis:
  \describe{
    \item{arclength}{The length of the test information or scale curve.}
    \item{arclengthvec}{Positions on the test information or scale curve corresponding to a
      fine mesh of score index values (typically 101 values between 0 and 100).} 
    \item{arclengthfd}{Functional data object representing the relation between
      the score index abscissa and the arclength or information ordinate.}
    \item{theta_al}{A vector of positions on the test information or scale curve corresponding
      to the input score index values in argument \code{theta}.} 
    \item{Qvec_al}{Values in arc length of the five marker percentages.}
    \item{binctr_al}{Values in arc length of the bin centers.}
    \item{Wfd.info}{A functional data object representing the relation between
      the arclength or information abscissa and the score index ordinate.}
    \item{Wdim.index}{The dimension of the overspace, which equal to sum of the number of 
      options in the items specified in \code{itemindex}.}
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
  \code{\link{Wpca.plot}}
}
\examples{
#  Example 1.  Display the arc length curve for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
indfine       <- seq(0,100,len=101)
WfdList       <- Quant_13B_problem_parList$WfdList
theta         <- Quant_13B_problem_parList$theta
Qvec          <- Quant_13B_problem_parList$Qvec
binctr        <- Quant_13B_problem_parList$binctr
arclenList    <- theta2arclen(theta, Qvec, WfdList, binctr)
names(arclenList)
oldpar <- par(no.readonly=TRUE)
ArcLength.plot(arclenList$arclength, arclenList$arclengthvec)
#  Example 2.  Display the arc length curve for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
par(oldpar)
}
