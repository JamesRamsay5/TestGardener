\name{entropies}
\alias{entropies}
\title{
Entropy measures of inter-item dependency
}
\description{
Entropy \eqn{I_1} is a scalar measure of how much information is required to predict 
the outcome of a choice number 1 exactly, and consequently is a measure of item effectiveness suitable for multiple choice tests and rating scales.  
Joint entropy \eqn{J_{1,2}} is a scalar measure of the cross-product of multinomial 
vectors 1 and 2. Mutual entropy \eqn{I_{1,2} = I_1 + I_2 - J_{1,2}} is a measure 
of the co-dependency of items 1 and 2, and thus the analogue of the negative
log of a squared correlation \eqn{R^2}.  this function computes all four types
of entropies for two specificed items.
}
\usage{
entropies(index, m, n, chcemat, noption)
}
\arguments{
  \item{index}{A vector of length N containing score index values for each
  test taker.}
  \item{m}{The index of the first choice.}
  \item{n}{The index of the second choice.}
  \item{chcemat}{The data matrix containing the indices of choisen options for
  each test taker.}
  \item{noption}{A vector containing the number of options for all items.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
  \item{I_m: }{The entropy of item m.}
  \item{I_n: }{The entropy of item n.}
  \item{J_nm: }{The joint entropy of items m and n.}
  \item{I_nm: }{The mutual entropy of items m and n.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Entropy_plot}}
}
\examples{
#  Load needed objects
chcemat <- Quant_13B_problem_dataList$chcemat
index   <- Quant_13B_problem_parmList$index
noption <- matrix(5,24,1)
#  compute mutual entropies for all pairs of the first 6 items
Mvec    <- 1:6
Mlen    <- length(Mvec)
Hmutual <- matrix(0,Mlen,Mlen)
for (i1 in 1:Mlen) {
  for (i2 in 1:i1) {
    Result <- entropies(index, Mvec[i1], Mvec[i2], chcemat, noption)
    Hmutual[i1,i2] = Result$Hmutual
    Hmutual[i2,i1] = Result$Hmutual
  }
}
print("Matrix of mutual entries (off-digagonal) and self-entropies (diagonal)")
print(round(Hmutual,3))
}