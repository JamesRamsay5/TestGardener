\name{Sensitivity_plot}
\alias{Sensitivity_plot}
\title{
  Plots all the sensitivity curves for selected items or questions.}
\description{
  A sensitivity curve for an option is the first derivative of the corresponding
  surprisal curve.  Its values can be positive or negative, and the size of the departure
  from zero at any point on the curve is the amount information contributed by that curve 
  to locating the value of an examinee or respondent on the score index continuum.
}
\usage{
Sensitivity_plot(scrfine, SfdList, Qvec, dataList, plotindex=1:n, 
                 plotrange=c(min(scrfine),max(scrfine)),   
                 key=NULL, titlestr=NULL, saveplot=FALSE, width=c(-0.2,0.2), 
                 ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
}
\arguments{
  \item{scrfine}{A vector of length nfine (usually 101) containing equally spaced
  	points spanning the plotrange.  Used for plotting.}
  \item{SfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  
    Its length is equal to the number of items in the test or questions in 
    the scale.  Each member of \code{SfdList} is a named list containing 
    information computed during the analysis.}
  \item{Qvec}{The values of the five marker percentiles.}
  \item{dataList}{A list that contains the objects needed to analyse the test 
    or rating scale.}
  \item{plotindex}{A set of integers specifying the numbers of the items or 
    questions to be displayed.}
  \item{plotrange}{A vector of length 2 containing the plot boundaries within or
    over the score index interval c(0,100).}
  \item{key}{A integer vector of indices of right answers.  
    If the data are rating scales, this can be NULL.}
  \item{titlestr}{A title string for plots.}
  \item{saveplot}{A logical value indicating whether the plot should be saved to a pdf file.}
  \item{width}{A vector of length 2 defining the lower and upper limits on the ordinate 
    for the plots.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
}
\details{
  Sensitivity curves for each question indexed in the \code{index} argument.  
  A request for a keystroke is made for each question.
}
\value{
  A list vector is returned which is of the length of argument \code{plotindex}.  
  Each member of the vector is a \code{gg} or \code{ggplot} object for the associated
  \code{plotindex value}.  Each plot can be displayed using the \code{print} command.
  The plots of item power are produced as a side value even if no output object is
  specified in the call to the function.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Power_plot}},
  \code{\link{Entropy_plot}},
  \code{\link{ICC_plot}}
}
\examples{
#  Example 1.  Display the option sensitivity curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees.
dataList <- Quant_13B_problem_dataList
SfdList  <- Quant_13B_problem_parmList$SfdList
Qvec     <- Quant_13B_problem_parmList$Qvec
scrfine  <- seq(0,100,len=101)
oldpar   <- par(no.readonly=TRUE)
Sensitivity_plot(scrfine, SfdList, Qvec, dataList, plotindex=1)
par(oldpar)
}