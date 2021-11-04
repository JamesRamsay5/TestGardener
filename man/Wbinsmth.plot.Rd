\name{Wbinsmth.plot}
\alias{Wbinsmth.plot}
\title{Plot probability and surprisal curves for a selection of test or scale 
       items.}
\usage{
  Wbinsmth.plot(binctr, Qvec, WfdList, dataList,
                twoplot=TRUE, ptsplot=TRUE, alltype=TRUE, 
                landscape=FALSE, saveplot=FALSE, plotindex=1:n, Wrng=c(0,5),
                ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
}
\description{
Wbinsmth.plots plots each item in argument \code{plotindex} in turn after 
function\code{Wbinsmth()} has used spline smoothing to estimate item and 
option characteristic curves.
}
\arguments{
  \item{binctr}{A vector of bin center values.}
  \item{Qvec}{A vector of five marker percentile values.}
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
      \item{chartList:}{A list vector of length \code{n}.  Each list contains the
      \code{K} by \code{M-1} matrix of initial values for the surprisal curves for
      the corresponding item.}
    }
  }
  \item{dataList}{A list that contains the objects needed to analyse the test or 
  rating scale with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
             N is the number of examinees or respondents and
             n is the number of items.
             Entries in the matrices are the indices of the options chosen.
             Column i of U is expected to contain only the integers 
             \code{1,...,noption}.}  
      \item{optList:}{A list vector containing the numerical score values 
      assigned to the options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is 
       	a numeric vector of length n containing the indices the right answers.  
       	Otherwise, it is NULL.}
      \item{chartList:}{A list vector length n each member of which is a 
      matrix of parameters values defining a set of M surprisal curves.  
      The matrix has K rows and M - 1 columns where K is the number of basis 
      functions defining the curves and M is the number curves.}
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
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing 
      the bin boundaries alternating with the bin centers. These are initially 
      defined as \code{seq(0,100,len=2*nbin+1)}.}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 
      5, 25, 50, 75 and 95.}
    }
  }
  \item{twoplot}{logical, if true plot both P and W curves.}
  \item{ptsplot}{logical, if true plot observed proportions.}
  \item{alltype}{logical, if true, plot right and wrong P curves.}
  \item{landscape}{logical, if true plots are side by side.}
  \item{plotindex}{A vector of indices of items to be plotted.}
  \item{saveplot}{A logical value indicating whether the plot should be saved to a pdf file.}
  \item{Wrng}{A vector of length 2 specifying the range for plotting surprisal values.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
}
\value{
A list vector is returned which is of the length of argument \code{plotindex}.  Each
member of the vector is a \code{gg} or \code{ggplot} object for the associated
\code{plotindex value}.  Each plot can be displayed using the \code{print} command.
The plots of item power are produced as a side value even if no output object is
specified in the call to the function.
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
  \code{\link{Sensitivity.plot}},
  \code{\link{Power.plot}},
  \code{\link{Entropy.plot}},
  \code{\link{Item.plot}},
  \code{\link{Wbinsmth}},
}
\examples{
#  Example 1.  Display the item surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList  <- Quant_parList$WfdList
Qvec     <- Quant_parList$Qvec
binctr   <- Quant_parList$binctr
titlestr <- NULL
#  plot the curves for the first question
Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, plotindex=1)
#  Example 2.  Display the item proability and surprisal curves for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}


