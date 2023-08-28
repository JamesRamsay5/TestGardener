\name{DHfun}
\alias{DHfun}
\title{Compute the first and second derivatives of the negative log likelihoods}
\description{
DHfun computes the first and second derivatives of the negative log likelihoods for 
a set of examinees.
Items can be either binary or multi-option.
The analysis is within the closed interval [0,100].
}
\usage{
  DHfun(theta, WfdList, Umat)
}
\arguments{
  \item{theta}{Initial values for score indices in [0,n]/[0,100]. Vector of size N.}
  \item{WfdList}{A numbered list object produced by a TestGardener analysis of a test.  Its length
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
  \item{Umat}{An \code{N} by \code{n} matrix of responses.  If \code{N} = 1, it 
  can be a vector of length \code{n}.}
}
\value{
A named list for results \code{DH} and \code{D2H}:
  \describe{
    \item{DH:}{ First derivatives of the negative log likelihood values, 
      vector of size N}
    \item{D2H:}{Second derivatives of the negative log likelihood values, 
      vector of size N}
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
  \code{\link{make.dataList},
        \link{Hfun},
        \link{Hfuns.plot}
  }
}
\examples{
#  Example 1:
#  Compute the first and second derivative values of the objective function for 
#  locating each examinee for the 24-item short form of the SweSAT quantitative  
#  test on thepercentile score index continuum.
WfdList <- Quant_13B_problem_parList$WfdList
theta   <- Quant_13B_problem_parList$theta
U       <- Quant_13B_problem_dataList$U
DHfunResult <- DHfun(theta, WfdList, U)
DHval  <- DHfunResult$DH
D2Hval <- DHfunResult$D2H
print(paste("Mean               of objective gradient =",round(mean(DHval),4)))
print(paste("Standard deviation of objective gradient =",round(sqrt(var(DHval)),4)))
print(paste("Mean               of objective Hessian =",round(mean(D2Hval),4)))
print(paste("Standard deviation of objective Hessian =",round(sqrt(var(D2Hval)),6)))
#  Example 2:
#  Compute the first and second derivative values of the objective function for 
#  locating each examinee for the 13-item Symptom Distress scale  
#  on the percentile score index continuum.
#  Proceed as above changing "Quant" for "SDS".
}
