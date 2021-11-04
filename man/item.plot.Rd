\name{Item.plot}
\alias{Item.plot}
\title{
  Four plots for item assessment
}
\description{
Plots of probability curves, sensitivity curves, power curves and
entropy curve for selected items.
}
\usage{
  Item.plot(WfdList, Qvec, dataList, plotindex=1:n, key=NULL,
            titlestr=NULL, saveplt=FALSE,
            ttlsz=NULL,axisttl=NULL,axistxt=NULL,lgdlab=NULL,
            width=c(-0.2,0.2), height=1)
}
\arguments{
  \item{WfdList}{List vector of length code{n} containing named list of
  objects specifying probability and surprisal curves.}
  \item{Qvec}{Score index values for current marker percentages.}
  \item{dataList}{Named list of objects required for analysis of data.}
  \item{plotindex}{Indices of item plots to be displayed.}
  \item{key}{Vector of length \code{n} of right answer option indices.  
  Only required for multiple choice test data.}
  \item{titlestr}{String for title of data.}
  \item{saveplt}{Logical value.  If TRUE, plots are saved in a pdf file.}
  \item{ttlsz}{Integer specifying font size for title.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
  \item{width}{Vector of length 2 of ordinate axis limits.}
  \item{height}{Value for shifting marker value labels.}
}
\details{
The four panels in these plots provide graphical information about the 
performance of items as indicators of optimal score index values.

The probability curves are of probabilities of choice functions of score index.

The sensitivity curves are of the derivatives of surprisal curves, and their
departure from 0 is a direct indication of the strength of the information
about whether a specific score index value should be reduced or increased.

The power curve is the root-sum-square of the sensitivity curves, and
indicate how much information at any score index value the item provides
for indicating the value of the score index value.

The entropy curve displays the level of randomness in the surprisal curves
at a score index value, which is maximized when all surprisal values 
are the same, and is minimized at zero when only one choice is non-zero.
For multiple choice tests, the entropy curve should decrease from one
on the right ot near zero on the left.
}
\value{
The plots are a side effect.  If \code{saveplt=TRUE, a pdf file is returned
containing all the plots for the selected items.}}
\references{
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
  Journal of Educational and Behavioral Statistics, 45, 297-315.

  Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
  information-based psychometrics.  Psych, 2, 347-360.

  http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Wbinsmth.plot}},
  \code{\link{Sensitivity.plot}},
  \code{\link{Power.plot}},
  \code{\link{Entropy.plot}}
}
