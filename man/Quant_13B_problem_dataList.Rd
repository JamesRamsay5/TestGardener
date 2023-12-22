\name{Quant_13B_problem_dataList}
\alias{Quant_13B_problem_dataList}
\docType{data}
\title{
  List of objects essential for an analysis of the abbreviated SweSAT 
  Quantitative multiple choice test.
}
\description{
  The data are for 1000 randomly selected examinees taking 24 math analysis 
  multiple choice questions.
}
\usage{Quant_13B_problem_dataList}
\format{A named list.}
\details{
  A named list with 19 members:
  \describe{
      \item{chcemat:}{A matrix of response data with N rows and n columns where
        N is the number of examinees or respondents and n is the number of items.
        Entries in the matrices are the indices of the options chosen.
        Column i of chcemat is expected to contain only the integers  
        \code{1,...,noption}.}  
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is 
       	a numeric vector of length n containing the indices the right answers.  
       	Otherwise, it is NULL.}
      \item{titlestr:}{A string containing a title for the analysis.}
      \item{N:}{The number of persons tested}
      \item{n:}{The number of questons or items}
      \item{noption:}{A numeric vector of length n containing the numbers of 
        options for each item.}
      \item{Sdim:}{The total dimension of the surprisal scores.}
      \item{grbgvec:}{A vector of length indicating which option for each item 
        contains missing or illegal choice values.  If 0, there is no such 
        option.}
      \item{ScoreList:}{A list vector or length n with each object a numeric 
        vector of weights assigned to each option for each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{NumBasis:}{The number of spline basis functions.}
      \item{Sbasis:}{An basisfd object for the defining the surprisal curves.}
      \item{itemlabvec:}{A character vector with a title string for each item.}
      \item{optlabList:}{A list vector of length n with a character vector of 
        labels for each object within each item.}
      \item{scrvec:}{A vector of length N containing the examinee or 
        respondent sum scores.}
      \item{itmvec:}{A vector of length n containing the item sum scores.}
      \item{scrjit:}{A numeric vector of length N containing small jitters to
        each sum score to break up ties,}
      \item{sumscr_rng:}{A vector of length 2 containing the limits of observed 
        sum scores.}
      \item{SfdList:}{A list vector containing essential objects for each item.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{indexQnt:}{A numeric vector of length 2*nbin + 1 containing 
        the bin boundaries alternating with the bin centers. These are initially 
        defined as \code{seq(0,100,len=2*nbin+1)}.}
      \item{percntrnk:}{A vector length N containing the sum score 
        percentile ranks.}
      \item{PcntMarkers:}{The marker percentages for plotting: 
        5, 25, 50, 75 and 95.}
    }
}
\keyword{datasets}
