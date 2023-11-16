\name{Quant_13B_problem_parmList}
\alias{Quant_13B_problem_parmList}
\docType{data}
\title{
  Parameter list for 24 items from the quantitative SweSAT subtest.
}
\description{
  The data are for 1000 examinees randomly selected from those who took the 2013 
  quantitative subtest of the SweSAT university entrance exam. 
  The questions are only the 24 math analysis questions, and each question has 
  four options.
  The analysis results are after 10 cycles of alternating between estimating 
  surprisal curves and estimating percentile score index values.
}
\usage{Quant_13B_problem_parmList}
\format{A named list.}
\value{
  The object \code{Quant_13B_problem_parmList} is a named list with these members:
  \item{index:}{A vector of length \code{N} of estimated values of the percentile 
    rank score index.}
  \item{indexQnt:}{A vector of length 2*nbin + 1 containing bin boundaries 
  	alternating with bin centres.}
  \item{SfdList:}{A list vector of length equal to the number of questions.  
    Each member contains eight results for the surprisal curves 
    associated with a question.}
  \item{logdensfd:}{A functional data object representing the logarithm of the 
    density of the percentile rank score index values.}
  \item{C:}{The norming constant: the density function is \code{exp(logdensfd)/C}.}
  \item{densfine:}{A fine mesh of probability density values of the percentile 
    rank score index.}
  \item{denscdf:}{A fine mesh of cumulative probability distribution values 
    used for interpolating values.}
  \item{Qvec:}{The score index values associated with the five marker percentages 
  	5, 25, 50, 75 and 95.}      
  \item{binctr:}{A vector of length nbin containing the centres of the bins.}
  \item{bdry:}{A vector of length nbin+1 containing the boundaries of the bins.}
  \item{freq:}{An nbin by M matrix of frequencies with which options are chosen.}
  \item{Smax:}{A maximum surprisal value used for plotting purposes.}
  \item{Hval:}{The value of the fitting criterion \code{H} for a single 
    examinee or respondent.}
  \item{DHval:}{The value of the first derivative of the fitting criterion 
    \code{H} for a single examinee or respondent.}
  \item{D2Hval:}{The value of the second derivative of the fitting criterion 
    \code{H} for a single examinee or respondent.}
  \item{active:}{A logical vector of length N indicating which estimates of 
  	index are converged (FALSE) or not converged (TRUE).}
  \item{infoSurp:}{The length in bits of the test information curve.}
  \item{infofine:}{A mesh of 101 equally spaced positions along the test 
  	information curve.}
  \item{Qinfovec:}{The positions of the five marker percentages on the test 
  	information curve.}
  \item{scopevec:}{A vector of length N containing the positions of each 
   	examinee or respondent on the test information curve.}
}
\keyword{datasets}
