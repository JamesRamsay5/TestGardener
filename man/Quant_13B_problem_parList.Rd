\name{Quant_13B_problem_parList}
\alias{Quant_13B_problem_parList}
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
\usage{Quant_13B_problem_parList}
\format{A named list.}
\value{
  The object \code{Quant_13B_problem_parList} is a named list with these members:
  \describe{
    \item{theta:}{A vector of length N of estimated values of the percentile 
    	rank score index.}
    \item{thetaQnt:}{A vector of length 2*nbin + 1 containing bin boundaries 
    	alternating with bin centres.}
    \item{WfdList:}{A list vector of length equal to the number of questions.  
        Each member contains the following eight results for the surprisal curves 
        associated with a question:
      \describe{
        \item{Wfd:}{A functional data object representing the M surprisal curves defined 
        	over the percentile rank range, where M is the number of options for that 
        	question.}
        \item{M:}{The number of option choices for the question.}
        \item{type:}{The type of question: 1 for multiple choice, 0 for rating scale.}
        \item{Pbin:}{A nbin by M matrix of proportions of examinees choosing each option 
        	for each bin,.}
        \item{Wbin:}{A nbin by M matrix of surprisal values of examinees choosing 
        	each option for each bin.}
        \item{Pmatfine:}{A nfine by M matrix of proportions over a fine mesh of 
        	proportions of examinees choosing each option.}
        \item{Wmatfine:}{A nfine by M matrix of surprisal values over a fine mesh 
        	of values of examinees choosing each option.}
        \item{DWmatfine:}{A nfine by M matrix of the first derivative of surprisal 
        	values over a fine mesh of values of examinees choosing each option.  
        	These are referred to as the sensitivity values.}
      }
    }
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
    \item{Wmax:}{A maximum surprisal value used for plotting purposes.}
    \item{Hval:}{The value of the fitting criterion \code{H} for a single 
    	examinee or respondent.}
    \item{DHval:}{The value of the first derivative of the fitting criterion 
    	\code{H} for a single examinee or respondent.}
    \item{D2Hval:}{The value of the second derivative of the fitting criterion 
    	\code{H} for a single examinee or respondent.}
    \item{active:}{A logical vector of length N indicating which estimates of 
    	theta are converged (FALSE) or not converged (TRUE).}
    \item{arclength:}{The length in bits of the test information curve.}
    \item{alfine:}{A mesh of 101 equally spaced positions along the test 
    	information curve.}
    \item{Qvec_al:}{The positions of the five marker percentages on the test 
    	information curve.}
    \item{theta_al:}{A vector of length N containing the positions of each 
    	examinee or respondent on the test information curve.}
  }
}
\keyword{datasets}
