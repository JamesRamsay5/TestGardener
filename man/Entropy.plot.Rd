\name{Entropy.plot}
\alias{Entropy.plot}
\title{
Plot item entropy curves for selected items or questions.
}
\description{
Item the value of the entropy curve at a point \code{theta} is the expected value of the surprisal curve values.  Entropy is a measure of the randomness of the surprisal value, which is
maximized when all the surprisal curves have the same value and has a minimum of zero if all
but a single curve has probability zero.  This is unattainable in the calculation, but can be
arbitrarily close to this state.}  
\usage{
  Entropy.plot(WfdList, Qvec, dataList, plotindex=1:n, height=1.0, value=0,  
               saveplot=FALSE, ttlsz=NULL, axisttl=NULL, axistxt=NULL)
}
\arguments{
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
  \item{Qvec}{The five marker percentile values.}
  \item{dataList}{A list vector containing objects essential to an analysis.}
  \item{plotindex}{
  A set of integers specifying the numbers of the items or questions to be displayed.
  }
  \item{height}{
  A positive real number defining the upper limit on the ordinate for the plots.
  }
  \item{value}{Number required by ggplot2.  Defaults to 0.}
  \item{saveplot}{A logical value indicating whether the plot should be saved to a pdf file.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
}
\details{
  An entropy curve for each question indexed in the \code{index} argument.  
  A request for a keystroke is made for each question.  The answer to question
  strongly defines the optimal position of an estimated score index value where
  the curve is high value.  Values of entropy curves typically range over [0,1].
}
\value{
The plots of the entropy curves specified in \code{plotindex} are produced as a side
effect.
If \code{saveplot} is TRUE, the plots of item entropy curves specified in \code{plotindex} are bundled into a single postscript or .pdf file and the file name is defined by 
\code{paste(dataList$titlestr,i,'-entropy.pdf',sep="")}.  The file is then output
as a returned value.
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
  \code{\link{Item.plot}},
  \code{\link{Hfuns.plot}},
  \code{\link{Wbinsmth.plot}}
}
\examples{
#  Example 1.  Display the item entropy curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  plot the entropy curve for the first item
WfdList <- Quant_parList$WfdList
Qvec    <- Quant_parList$Qvec
Entropy.plot(WfdList, Qvec, Quant_dataList, plotindex=1)
#  Example 2.  Display the item entropy curves for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}