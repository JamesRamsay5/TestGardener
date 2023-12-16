\name{mu}
\alias{mu}
\title{Compute the expected test score by substituting probability of choices
for indicator variable 0-1 values.
Binary items assumed coded as two choice items.}
\usage{
  mu(index, SfdList, scoreList)
}
\arguments{
  \item{index}{Initial values for score indices in the interval [0,100]. 
    A vector of size \code{N}.}
  \item{SfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  Its length
    is equal to the number of items in the test or questions in the scale.  
    Each member of \code{SfdList} is a named list containing information computed 
    during the analysis.}
  \item{scoreList}{A numbered list of length n.  Each member contains the weights
    assigned to each option for that item or question.}
}
\description{
  Compute the expected test score by substituting probability of choices
  for indicator variable 0-1 values.
  Binary items assumed coded as two choice items.
}
\value{A vector of test score values.}
\author{Juan Li and James Ramsay}
\references{
Ramsay, J. O., Li J. and Siberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Siberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\seealso{
  \code{\link{scoreDensity}}
}
\examples{
#  Example 1.  Compute expected sum score values for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
scoreList <- Quant_13B_problem_dataList$scoreList
SfdList   <- Quant_13B_problem_parmList$SfdList
index     <- Quant_13B_problem_parmList$index
muvec     <- mu(index, SfdList, scoreList)
par(c(1,1))
hist(muvec,11)
}
