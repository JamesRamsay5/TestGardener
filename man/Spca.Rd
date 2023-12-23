\name{Spca}
\alias{Spca}
\title{
Functional principal components analysis of information curve
}
\description{
  A test or scale analysis produces a space curve that varies with in the space of 
  possible option curves of dimension \code{Sdim}.  Fortunately, it is usual that most
  of the shape variation in the curve is within only two or three dimensions, and these
  can be fixed by using functional principal components analysis.
}
\usage{
  Spca(SfdList, nharm=2, Sdim=NULL, rotate=TRUE)
}
\arguments{
  \item{SfdList}{A numbered list object produced by a TestGardener analysis 
  of a test.  Its length is equal to the number of items in the test or 
  questions in the scale.  Each member of \code{SfdList} is a named list 
  containing information computed during the analysis.}
  \item{Sdim}{Interval over which curve is plotted. All if Sdim == NULL.}
  \item{nharm}{The number of principal components of the test information or scale 
  curve to be used to display the curve.  Must be either 2 or 3.}
  \item{rotate}{If true, rotate principal components of the test information or 
  scale curve to be used to display the curve to VARIMAX orientation.}
}
\value{
A named list with these members:
  \item{harmvarmxfd}{Functional data objects for the principal
    components of the curve shape.} 
  \item{varpropvarmx}{Proportions of variance accounted for by the principal 
    components}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Spca_plot}}
}
\examples{
#  Example 1.  Display the test information curve for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  plot a two-dimension version of manifold curve
Sdim     <- Quant_13B_problem_dataList$Sdim
SfdList  <- Quant_13B_problem_parmList$SfdList
index    <- Quant_13B_problem_parmList$index
infoSurp <- Quant_13B_problem_parmList$infoSurp 
#      <- Quant_13B_problem_dataList$Sdim
oldpar <- par(no.readonly=TRUE)
on.exit(oldpar)
Results <- Spca(SfdList, nharm=2, rotate=FALSE)
varprop <- Results$varpropvarmx
print("Proportions of variance accounted for and their sum:")
print(round(100*c(varprop,sum(varprop)),1))
#  plot a three-dimension version of manifold curve
SfdList   <- Quant_13B_problem_parmList$SfdList
index     <- Quant_13B_problem_parmList$index
infoSurp  <- Quant_13B_problem_parmList$infoSurp 
Results   <- Spca(SfdList, nharm=3, rotate=FALSE)
varprop   <- Results$varpropvarmx
print("Proportions of variance accounted for and their sum:")
print(round(100*c(varprop,sum(varprop)),1))
}
