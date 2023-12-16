\name{DFfun}
\alias{DFfun}
\title{Compute the first and second derivatives of the negative log likelihoods}
\description{
DFfun computes the first and second derivatives of the negative log likelihoods for 
a set of examinees.
Items can be either binary or multi-option.
The analysis is within the closed interval [0,100].
}
\usage{
  DFfun(index, SfdList, chcemat)
}
\arguments{
  \item{index}{Initial values for score indices in [0,n]/[0,100]. 
  Vector of size N.}
  \item{SfdList}{A numbered list object produced by a TestGardener analysis of 
  a test.  Its length is equal to the number of items in the test or questions 
  in the scale. Each member of \code{SfdList} is a named list containing 
  information computed during the analysis.}
  \item{chcemat}{An \code{N} by \code{n} matrix of responses.  If \code{N} = 1, it 
  can be a vector of length \code{n}.}
}
\value{
  A named list for results \code{DF} and \code{D2F}:
  \item{DF:}{ First derivatives of the negative log likelihood values, 
   vector of size N}
  \item{D2F:}{Second derivatives of the negative log likelihood values, 
   vector of size N}
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
  #  Example 1:
  #  Compute the first and second derivative values of the objective function  
  #  for locating each examinee for the 24-item short form of the  
  #  SweSAT quantitative test on the percentile score index continuum.
  #  Use only the first five examinees.
  chcemat <- Quant_13B_problem_dataList$chcemat
  SfdList <- Quant_13B_problem_parmList$SfdList
  index   <- Quant_13B_problem_parmList$index
  DFfunResult <- DFfun(index[1:5], SfdList, chcemat[1:5,])
  DFval  <- DFfunResult$DF
  D2Fval <- DFfunResult$D2F
}
