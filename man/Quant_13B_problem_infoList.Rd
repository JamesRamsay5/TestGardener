\name{Quant_13B_problem_infoList}
\alias{Quant_13B_problem_infoList}
\docType{data}
\title{
  Arclength or information arameter list for 24 items from the quantitative 
  SweSAT subtest.
}
\description{
  The data are for 1000 examinees randomly selected from those who took the 2013 
  quantitative subtest of the SweSAT university entrance exam. 
  The questions are only the 24 math analysis questions, and each question has 
  four options.
  The analysis results are after 10 cycles of alternating between estimating 
  surprisal curves and estimating percentile score index values.
  The objects in list object \code{Quant_13B_problem_infoList} are required for plotting
  results over the arc length or information domain rather the score index domain.
  This domain is preferred because such plots are invariant with respect to changes 
  in the score index domain. It also has a metric structure so that differences 
  are comparable no matter where they fall within the information domain.  
}
\usage{Quant_13B_problem_infoList}
\format{A named list containing eight objects.}
\value{
  The object \code{Quant_13B_problem_parList} is a named list with these members:
  \describe{
    \item{arclength:}{The total length of the information domain measured in 
    	M-bits, where \code{M } is the number of options for a question.}
    \item{Wfd.theta:}{The log derivative functional data object defining 
    	a strictly increasing set of arc length values
        corresponding to set of score index values.}
    \item{arclengthvec:}{A mesh of equally-spaced values of 
    	indefinite integrals of sum of norms of surprisal derivatives.}
    \item{theta_al}{The N arc length values corresponding to the N
    	estimated score index values assigned to N examinees.}
    \item{Qvec_al:}{The arc length positions corresponding to the marker 
    	percentages 5, 25, 50, 75 and 95.}      
    \item{Wfd.info:}{The log derivative functional data object defining
        a strictly increasing set of score index valuescorresponding to 
        a set of arc length values.}
    \item{thetavec:}{A vector of score index values resulting from using 
    	function \code{monfd} with equally spaced arc length values and 
        \code{Wfd.info}.}
    \item{Wdim:}{The dimension of the over space containing the  surprisal 
    	pcurves.}
  }
}
\keyword{datasets}
