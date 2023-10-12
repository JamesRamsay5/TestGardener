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
  Wpca.plot(harmvarmxfd, nharm=2, titlestr=NULL)
}
\arguments{
  \item{harmvarmxfd}{Functional data objects for the principal
    components of the curve shape.} 
  \item{nharm}{Number of principal components.} 
  \item{titlestr}{A string for the title of the plot.  Defaults to NULL. }
}
\value{
Side effect is a two or three-dimensional plot of the principal component
approximation of the information curve using the plotly package.  Function
\code{plot_ly} does not return a value, but does render the graphic.
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
  \code{\link{Wpca}}
}
\examples{
#  Example 1.  Display the test information curve for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  plot a two-dimension version of manifold curve
WfdList    <- Quant_13B_problem_parList$WfdList
theta      <- Quant_13B_problem_parList$theta
arclength  <- Quant_13B_problem_parList$arclength 
Results    <- Wpca(WfdList,  nharm=2, rotate=TRUE)
varprop    <- Results$varpropvarmx
titlestr   <- "SweSAT problem items"
oldpar <- par(no.readonly=TRUE)
on.exit(oldpar)
Wpca.plot(Results$harmvarmxfd, nharm=2, titlestr)
print("Proportions of variance accounted for and their sum:")
print(round(100*c(varprop,sum(varprop)),1))
#  plot a three-dimension version of manifold curve
WfdList    <- Quant_13B_problem_parList$WfdList
theta      <- Quant_13B_problem_parList$theta
arclength  <- Quant_13B_problem_parList$arclength 
Results    <- Wpca(WfdList,  nharm=3, rotate=TRUE)
varprop    <- Results$varpropvarmx
Wpca.plot(Results$harmvarmxfd, nharm=3, titlestr)
print("Proportions of variance accounted for and their sum:")
print(round(100*c(varprop,sum(varprop)),1))
#  Example 2.  Display the test information curve for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}
