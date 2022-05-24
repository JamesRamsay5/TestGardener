\name{scorePerformance}
\alias{scorePerformance}
\title{
  Calculate mean squared error and bias for a set of score index values from 
  simulated data.
}
\description{
  This function is used on \code{ConditionalSimulation} after the simulated 
  data matrices have been analyzed.  It prepares the objects necessary for the 
  performance plots produced by functions \code{RMSEbias1.plot} and 
  \code{RMSEbias2.plot}.
}
\usage{
  scorePerformance(dataList, simList)
}
\arguments{
  \item{dataList}{A list that contains the objects needed to analyse the test 
  or rating scale with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
        N is the number of examinees or respondents and n is the number of items.
        Entries in the matrices are the indices of the options chosen.
        Column i of U is expected to contain only the integers 
         \code{1,...,noption}.}  
      \item{optList:}{A list vector containing the numerical score values 
        assigned to the options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is 
       	a numeric vector of length n containing the indices the right answers.  
       	Otherwise, it is NULL.}
      \item{WfdPar:}{An fdPar object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of 
        options for each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed 
        sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or 
        respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item 
        sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score 
        percentile ranks.}
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing the 
        bin boundaries alternating with the bin centers. These are initially 
        defined as \code{seq(0,100,len=2*nbin+1)}.}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 
        5, 25, 50, 75 and 95.}
    }
  }
  \item{simList}{A named list containing these objects:
    \describe{
      \item{sumscr:}{A matrix with row dimension \code{ntheta}, the number of 
        population score index values and column dimension \code{nsample}, the 
        number of simulated samples.}
      \item{theta:}{An \code{ntheta} by \code{nsample} of estimated score 
        index values.}
      \item{mu:}{An \code{ntheta} by \code{nsample} of estimated expected 
        score values.}
      \item{al:}{An \code{ntheta} by \code{nsample} of estimated test 
        information curve values.}
      \item{thepop:}{A vector of population score index values.}
      \item{mupop:}{A vector of expected scores computed from the population 
        score index values.}
      \item{alpop:}{A vector of test information values computed from the 
        population score index values.}
      \item{n:}{The number of questions.}
      \item{Qvec:}{The five marker percentile values.}
    }
  }
}
\value{
  A named list containing these objects:
    \describe{
      \item{sumscr:}{A matrix with row dimension \code{ntheta}, the number of 
        population score index values and column dimension \code{nsample}, the 
        number of simulated samples.}
      \item{theta:}{An \code{ntheta} by \code{nsample} matrix of estimated score 
        index values.}
      \item{mu:}{An \code{ntheta} by \code{nsample} matrix of estimated expected 
        score values.}
      \item{al:}{An \code{ntheta} by \code{nsample} matrix of estimated test 
        information curve values.}
      \item{thepop:}{A vector of population score index values.}
      \item{mupop:}{A vector of expected scores computed from the population 
        score index values.}
      \item{alpop:}{A vector of test information values computed from the 
        population score index values.}
      \item{n:}{The number of questions.}
      \item{Qvec:}{The five marker percentile values.}
    }
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
