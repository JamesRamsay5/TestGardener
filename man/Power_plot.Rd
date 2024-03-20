\name{Power_plot}
\alias{Power_plot}
\title{Plot item power curves for selected items or questions.}
\description{
  Item surprisal power curves are the square root of the sum over options of 
  the squared surprisal sensitivity curves.}
\usage{
  Power_plot(scrfine, SfdList, Qvec, dataList, plotindex=1:n, 
             plotrange=c(min(scrfine),max(scrfine)), height=0.5,  
             value=0, ttlsz=NULL, axisttl=NULL, axistxt=NULL)
}
\arguments{
  \item{scrfine}{A vector of length nfine (usually 101) containing equally spaced
  	points spanning the plotrange.  Used for plotting.}
  \item{SfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  
    Its length is equal to the number of items in the test or questions in 
    the scale.  Each member of \code{SfdList} is a named list containing 
    information computed during the analysis.}
  \item{Qvec}{The five marker percentile values.}
  \item{dataList}{A list vector containing objects essential to an analysis.}
  \item{plotindex}{A set of integers specifying the numbers of the items or 
  questions to be displayed.}
  \item{plotrange}{A vector of length 2 containing the plot boundaries within or
    over the score index interval c(0,100).}
  \item{height}{A positive real number defining the upper limit on the ordinate 
  	for the plots.}
  \item{value}{Number required by ggplot2.  Defaults to 0.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
}
\details{
  A surprisal power curve for each question indexed in the \code{index} argument.  
  A request for a keystroke is made for each question.  The answer to question
  strongly defines the optimal position of an estimated score index value where
  the curve is high value.  Values of power curves typically range over [0,0.5].
}
\value{
  The plots of the power curves specified in \code{plotindex} are produced as 
  a side effect.
  If \code{saveplot} is TRUE, the plots of item power curves specified in 
  \code{plotindex} are bundled into a single postscript or .pdf file and the 
  file name is defined by \code{paste(dataList$titlestr,i,'-power.pdf',sep="")}.  
  The file is then outputas a returned value.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Sensitivity_plot}},
  \code{\link{Entropy_plot}},
  \code{\link{Ffuns_plot}},
  \code{\link{ICC_plot}}
}
\examples{
#  Example 1.  Display the item power curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  plot the power curve for the first item
dataList <- Quant_13B_problem_dataList
SfdList  <- Quant_13B_problem_parmList$SfdList
Qvec     <- Quant_13B_problem_parmList$Qvec
scrfine <- seq(0,100,len=101)
oldpar <- par(no.readonly=TRUE)
Power_plot(scrfine, SfdList, Qvec, dataList, plotindex=1)
par(oldpar)
}