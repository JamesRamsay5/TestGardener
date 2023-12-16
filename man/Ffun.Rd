\name{Ffun}
\alias{Ffun}
\title{Compute the negative log likelihoods associated with a vector of score index values.}
\description{
Ffun computes the negative log likelihoods for a set of examinees, each 
at a single value index.
}
\usage{
  Ffun(index, SfdList, chcemat)
}
\arguments{
  \item{index}{A vector of size \code{N} containing values for score indices in 
  the interval [0,100].}
  \item{SfdList}{A numbered list object produced by a TestGardener analysis of 
  a test.  Its length is equal to the number of items in the test or questions 
  in the scale.  Each member of \code{SfdList} is a named list containing 
  information computed during the analysis.}
  \item{chcemat}{An \code{N} by \code{n} matrix of responses or, for a single examinee,
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
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{make_dataList},
        \link{index_fun},
        \link{Ffun},
        \link{Ffuns_plot}
  }
}
\examples{
  #  Example 1:  Compute the values of the objective function for locating each  
  #  examinee or respondent for the 24-item short form of the SweSAT quantitative  
  #  test on the percentile score index continuum [0,100].
  #  Use only the first five examinees
  chcemat <- Quant_13B_problem_dataList$chcemat
  SfdList <- Quant_13B_problem_parmList$SfdList
  index   <- Quant_13B_problem_parmList$index
  Fval    <- Ffun(index[1:5], SfdList, chcemat[1:5,])
}