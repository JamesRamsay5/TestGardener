\name{Wpca.plot}
\alias{Wpca.plot}
\title{
Plot the test information or scale curve in either two or three dimensions.
}
\description{
  A test or scale analysis produces a space curve that varies with in the space of 
  possible option curves of dimension \code{Wdim}.  Fortunately, it is usual that most
  of the shape variation in the curve is within only two or three dimensions, and these
  can be fixed by using functional principal components analysis.
}
\usage{
  Wpca.plot(arclength, WfdList, Wdim=NULL, nharm=2, rotate=TRUE, 
  titlestr=NULL)
}
\arguments{
  \item{arclength}{The total length of the test information or scale curve as computed by function 
  \code{theta2arclen.}}
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
  \item{Wdim}{The total number of options in the test or scale.}
  \item{nharm}{The number of principal components of the test information or scale 
  curve to be used to display the curve.  Must be either 2 or 3.}
  \item{rotate}{If true, rotate principal components of the test information or 
  scale curve to be used to display the curve to VARIMAX orientation.}
  \item{titlestr}{A string for the title of the plot.  Defaults to NULL. }
}
\value{
A named list with these members:
  \item{pcaplot}{If two dimensions or harmonics are specified, this is a 
  \code{gg} or \code{ggplot} object that can be displayed using the \code{print} 
  command.  If three dimensions are specified, this is NULL.}
  \item{harmvarmxfd}{Functional data objects for the principal
    components of the curve shape.} 
  \item{varpropvarmx}{Proportions of variance accounted for by the principal 
    components}
  The principal components are VARIMAX rotated by default.  The plot is displayed
  as a side value even if no output object is specified in the call to the function.
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
#  Example 1.  Display the test information curve for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  plot a two-dimension version of manifold curve
WfdList    <- Quant_13B_problem_parList$WfdList
theta      <- Quant_13B_problem_parList$theta
arclength  <- Quant_13B_problem_parList$arclength 
oldpar <- par(no.readonly=TRUE)
on.exit(oldpar)
Wpca.plotResults <- Wpca.plot(arclength, WfdList, Quant_13B_problem_dataList$Wdim)
varprop <- Wpca.plotResults$varpropvarmx
print("Proportions of variance accounted for and their sum:")
print(round(c(varprop,sum(varprop)),3))
#  plot a three-dimension version of manifold curve
WfdList    <- Quant_13B_problem_parList$WfdList
theta      <- Quant_13B_problem_parList$theta
arclength  <- Quant_13B_problem_parList$arclength 
Wpca.plotResults <- Wpca.plot(arclength, WfdList, Quant_13B_problem_dataList$Wdim, nharm=3)
varprop <- Wpca.plotResults$varpropvarmx
print("Proportions of variance accounted for and their sum:")
print(round(c(varprop,sum(varprop)),3))
#  Example 2.  Display the test information curve for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}
