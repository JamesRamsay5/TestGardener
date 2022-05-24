\name{Hfun}
\alias{Hfun}
\title{Compute the negative log likelihoods associated with a vector of score index values.}
\description{
Hfun computes the negative log likelihoods for a set of examinees, each 
at a single value theta.
}
\usage{
  Hfun(theta, WfdList, Umat)
}
\arguments{
  \item{theta}{A vector of size \code{N} containing values for score indices in the interval [0,100].}
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
  \item{Umat}{An \code{N} by \code{n} matrix of responses or, for a single examinee,
    a vector of length \code{n}.}
}
\value{
  A vector of length \code{N} of negative log likelihood values.
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
  \code{\link{make.dataList},
        \link{Hfun},
        \link{Hfuns.plot}
  }
}
\examples{
#  Example 1:  Compute the values of the objective function for locating each  
#  examinee or respondent for the 24-item short form of the SweSAT quantitative  
# test on the percentile score index continuum [0,100].
WfdList <- Quantshort_parList$WfdList
theta   <- Quantshort_parList$theta
U       <- Quantshort_dataList$U
Hval    <- Hfun(theta, WfdList, U)
print(paste("Mean               of objective function =",round(mean(Hval),2)))
print(paste("Standard deviation of objective function =",round(sqrt(var(Hval)),2)))
#  Example 2:  Compute the values of the objective function for locating each  
#  examinee or respondent for the 13-item Symptom Distress Scale  
# on the percentile score index continuum [0,100].
#  Proceed as above changing "Quant" for "SDS".
}