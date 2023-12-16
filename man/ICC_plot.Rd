\name{ICC_plot}
\alias{ICC_plot}
\title{Plot probability and surprisal curves for test or scale items.}
\usage{
  ICC_plot(scrfine, SfdList, dataList, Qvec, 
           binctr=NULL, data_point = FALSE, ci = FALSE, 
           plotType="S", Srng=c(0,5), DSrng=c(-0.2, 0.2), plotindex=1:n, 
           titlestr = NULL, itemscopevec = rep(0, length(plotindex)),  
           plotTitle = TRUE,autoplot = FALSE, plotMissing = TRUE, 
           plotrange=c(min(scrfine),max(scrfine)), shaderange = NULL,  
           ttlsz = NULL, axisttl = NULL, axistxt = NULL, 
           lgdlab = NULL, lgdpos = "bottom")
}
\description{
ICC_plot plots each item in argument \code{plotindex} in turn after 
function\code{Sbinsmth()} has used spline smoothing to estimate item and 
option characteristic curves.
}
\arguments{
  \item{scrfine}{A vector of 101 plotting points.}
  \item{SfdList}{A numbered list object produced by a TestGardener analysis of 
    a test.  Its length is equal to the number of items in the test or 
    questions in the scale.  Each member of \code{SfdList} is a named list 
    containing information computed during the analysis.}
  \item{dataList}{A list that contains the objects needed to analyse the test or 
    rating scale.}
  \item{Qvec}{A vector of five marker percentile values.  For
    plotting over information, this is replaced by Qinfovec returned
    as parmList$Qinfovec.}
  \item{binctr}{A vector of bin center values.  If the plot is over
    arc length or information, \code{binctr} is modified before calling
  \code{Sbinsth_plot} by the command
  \code{binctrinfo = pracma::interp1(indfine, alfine, binctr)}, and
    argument \code{binctr} is replaced by \code{binctrinfo}.}
  \item{data_point}{A logical value indicating whether to plot the data points.}
  \item{ci}{A logical value indicating whether to plot the confidence limits.}
  \item{plotType}{Type(s) of plot, default as "P" for probability, 
    can also be "S" for surprisal, 
    "DS" for sensitivity, and any combination of the three}
  \item{Srng}{A vector of length 2 specifing the plotting range for surprisal 
    values.}
  \item{DSrng}{A vector of length 2 specifing the plotting range for 
    sensitivity values.}
  \item{plotindex}{A vector of indices of items to be plotted.}
  \item{titlestr}{plot title}
  \item{itemscopevec}{A numeric vector containing item scope values.}
  \item{plotTitle}{indicator of showing the plot title, default as TRUE}
  \item{autoplot}{indicator for plotting all items in a batch}
  \item{plotMissing}{Determine if plot the extra option for missing/spoiled responses.}
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
  A list vector is returned which is of the length of argument \code{plotindex}.  
  Each member of the vector is a \code{gg} or \code{ggplot} object for the 
  associated \code{plotindex value}.  Each plot can be displayed using the 
  \code{print} command.
  The plots of item power are produced as a side value even if no output object 
  is specified in the call to the function.
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
  \code{\link{Power_plot}},
  \code{\link{Entropy_plot}},
  \code{\link{Sbinsmth}}
}
\examples{
#  Example 1.  Display the item surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
dataList <- Quant_13B_problem_dataList
SfdList  <- Quant_13B_problem_parmList$SfdList
Qvec     <- Quant_13B_problem_parmList$Qvec
binctr   <- Quant_13B_problem_parmList$binctr
infoSurpvec  <- Quant_13B_problem_infoList$infoSurpvec
Qinfovec     <- Quant_13B_problem_infoList$Qinfovec
bininfoctr   <- Quant_13B_problem_infoList$bininfoctr
titlestr <- "Quant_13B_problem"
#  plot the curves for the first question over the score index
oldpar  <- par(no.readonly=TRUE)
indfine <- seq(0,100,len=101)
ICC_plot(indfine, SfdList, dataList, Qvec, binctr, 
         data_point = TRUE, plotType = c("S", "P"), 
         Srng=c(0,4), plotindex=1)
#  plot the curves for the first question over test information
ICC_plot(infoSurpvec, SfdList, dataList, Qinfovec, bininfoctr,
          data_point = TRUE, plotType = c("S", "P"), 
          Srng=c(0,4), plotindex=1)
par(oldpar)
}


