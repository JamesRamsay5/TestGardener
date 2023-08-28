\name{plotICC}
\alias{plotICC}
\title{The most essential part for plotting probability/sueprisal
       curves of a single item.}
\usage{
  plotICC(Mi, scrfine, fitfinei, Qvec, keyi, 
          binctr, bin1, StdErr, 
          range, intercept, 
          ttllab, xlab, ylab, optVec,
          plotrange, shaderange, 
          ttlsz, axisttl, axistxt,
          lgdlab, lgdpos)
}
\description{
plotICC is the most essential part for plotting probability/sueprisal
curves of a single item..
}
\arguments{
  \item{Mi}{Number of options.}
  \item{scrfine}{A vector of 101 plotting points.}
  \item{fitfinei}{A matrix of probabilities/surprisal/sensitivity over a fine mesh.}
  \item{Qvec}{A vector of five marker percentile values.  For
    plotting over information, this is replaced by Qvec_al returned
    as parList$Qvec_al.}
  \item{keyi}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0.  
       	Otherwise, it is NULL.}
  \item{binctr}{A vector of bin center values.  If the plot is over
    arc length or information, \code{binctr} is modified before calling
  \code{Wbinsth.plot} by the command
  \code{binctr_al = pracma::interp1(indfine, alfine, binctr)}, and
    argument \code{binctr} is replaced by \code{binctr_al}.}
  \item{bin1}{A matrix containing proportions or surprisal values at each bin.}
  \item{StdErr}{A matrix of the values of the standard error of 
        probability/surprisal curves at each bin.}
  \item{range}{A vector of length 2 specifing the plotting range for y-axis.}
  \item{intercept}{A value indicating the y-intercept of the dashed horizontal line.}
  \item{ttllab}{plot title}
  \item{xlab}{title of x-axis, default as "Score Index"}
  \item{ylab}{title of y-axis, default as "Probability", Note: should work with plotType to define ylab for each subplot}
  \item{optVec}{vector of option labels for the legend}
  \item{plotrange}{A vector of length 2 containing the plot boundaries of
  	 the score index interval.}
  \item{shaderange}{a list of length 2 vector(s); set if users want to gray 
    out specific score range(s)}      
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
  \item{lgdpos}{legend position, could be set as "None" to remove the legend.}
}
\value{
  A ggplot object.
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
  \code{\link{Wbinsmth}},
}
\examples{
#  Example 1.  Display the item surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList  <- Quant_13B_problem_parList$WfdList
Qvec     <- Quant_13B_problem_parList$Qvec
binctr   <- Quant_13B_problem_parList$binctr
titlestr <- NULL
#  plot the curves for the first question
scrfine <- seq(0,100,len=101)
oldpar <- par(no.readonly=TRUE)
ICC.plot(scrfine, WfdList, Quant_13B_problem_dataList, Qvec, binctr, plotindex=1:2)
par(oldpar)
#  Example 2.  Display the item proability and surprisal curves for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}


