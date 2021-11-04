\name{SDS_dataList}
\alias{SDS_dataList}
\docType{data}
\title{
  List of objects essential for an analysis of the Symptom Distress Scale.
}
\description{
  The data are for 473 respondents to the Symptom Distress Scale.  Each question has
  rating scores 0, 1, 2, 3, and 4.
}
\usage{data("SDS_dataList")}
\format{A named list.}
\details{
  A named list with 15 members. See below for a description of each member.
  \describe{
    \item{U:}{An N by n matrix of integers indicating which option was chosen by one of 
    N examinees for a question.  The object should be in the matrix class and the 
    contents should be integers.  
    For any question the integers should be within the range 1 to the number of options.}
    \item{optList:}{A list vector object of length equal to the number of questions.  
    Each member of the list is a numeric vector of score values assigned to each option 
    by the test designer.  For multiple choice questions, these scores are usually 0 for 
    wrong answers and 1 for the right answer.  For rating scales they are often a 
    sequence of signed integers indicating intensity of some experience.}
    \item{grbg:}{A vector of length indicating which option for each item contains 
    missing or illegal choice values.  If 0, there is no such option.}
    \item{WfdPar:}{A functional parameter object that is required to represent the 
    surprisal curves. The number of basis functions should be small for smaller sizes and
    larger for large samples.  For example, 7 for N <= 500, 11 for N in [500,2000], 24 
    for N >= 2000.}
    \item{logdensfd:}{A functional data object used to represent the log density of the 
    percentile score index values.  The number of basis functions is often the same as 
    used for the WfdPar object.}
    \item{noption:}{An integer vector of length equal to the number of questions 
    indicating the number of options for each question.}
    \item{nbin:}{The number of bins used to contain the score index values.  The number 
    should ensure that each bin can hold at least 25 score values.  For medium sized 
    samples, 50 per bin is fine, and for large samples as many as are required to define 
    the detail in the surprisal curves.}
    \item{key:}{An integer vector of length equal to the number of questions and 
    containing the right answer indices.  This is not used for rating scales and can be 
    NULL in that case.}
    \item{scrrng:}{A vector of length 2 defining a range over which the sum and expected
    scores will be estimated.}
    \item{scrvec:}{A vector of length N containing the sum scores for the examinees or 
    respondents.}
    \item{percntrnk:}{The percentile ranks for each of the examinees or respondents.  
    This is the percentage of examinees with sum scores below that for an examinee.  Ties
    are often broken by adding a small random quantity to each sum score before computing
    ranks.}
    \item{scrfine:}{A fine mesh of score values used for plotting purposes, and running
    from scrrng[1] to scrrng[2].  Common number for mesh values are 51 or 101.}
    \item{scrtype:}{If the question is of multiple choice, the value is 1, and if not, 0
    .}
    \item{indrng:}{Defines the range of the score index.  Usually this is c(0,100).}
    \item{indfine:}{A fine mesh of score index values used for plotting, and usually 
    ranging from 0 to 100.  The number of mesh values is commonly 51 or 101.}
  }
}
\keyword{datasets}
