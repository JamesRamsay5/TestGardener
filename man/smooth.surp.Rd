\name{smooth.surp}
\alias{smooth.surp}
\title{
  Fit data with surprisal smoothing.
}
\description{
  Surprisal is -log(probability) where the logarithm is to the base being the 
  dimension \code{M} of the multinomial observation vector.  The surprisal 
  curves for each question are estimated by fitting the surprisal values of 
  binned data using curves whose values are within the \code{M-1}-dimensional 
  surprisal subspace that is within the space of non-negative M-dimensional 
  vectors.
}
\usage{
smooth.surp(argvals, y, Bmat0, Sfd, Zmat, wtvec=NULL, conv=1e-4,
            iterlim=50, dbglev=0)
}
\arguments{
  \item{argvals}{
    Argument value array of length N, where N is the number of observed
    curve values for each curve.  It is assumed that that these argument
    values are common to all observed curves.  If this is not the case,
    you will need to run this function inside one or more loops,
    smoothing each curve separately.
  }
  \item{y}{
    A \code{nbin} by \code{M_i} matrix of surprisal values to be fit.  
  }
  \item{Bmat0}{A \code{Snbasis} by \code{M_i - 1} matrix containing starting
  values for the iterative optimization of the least squares fit of the 
  surprisal curves to the surprisal data.
  }
  \item{Sfd}{
    A functional data object used toinitialize the optimization process.  
  }
  \item{Zmat}{An \code{M by M-1} matrix satisfying \code{Zmat'Zmat <- I} and 
    \code{Zmat'1 <- 0}.}
  \item{wtvec}{
    A vector of weights to be used in the smoothing.
  }
  \item{conv}{
    A convergence criterion.
  }
  \item{iterlim}{
    the maximum number of iterations allowed in the minimization of
    error sum of squares.
  }
  \item{dbglev}{
    Either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  If either level 1 or 2 is specified, it
    can be helpful to turn off the output buffering feature of S-PLUS.
  }
}
\value{
  A named list of class \code{surpFd} with these members:
  \item{PENSSE}{The final value of the penalized fitting criterion.}
  \item{DPENSSE}{The final gradient of the penalized fitting criterion.}
  \item{D2PENSSE}{The final hessian of the fitting criterion.}
  \item{SSE}{The final value of the error sum of squares.}
  \item{DSSE}{The final gradient of the error sum of squares.}
  \item{D2SSE}{The final hessian of the error sum of squares.}
  \item{DvecSmatDvecB}{The final cross derivative DvecSmatDvecX times 
    DvecXmatDvecB of the surprisal curve and the basis coordinates.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{eval.surp}},
  \code{\link{ICC_plot}},
  \code{\link{Sbinsmth}}
}
\examples{
  oldpar <- par(no.readonly=TRUE)
  #  evaluation points
  x <- seq(-2,2,len=11)
  #  evaluate a standard normal distribution function
  p <- pnorm(x)
  #  combine with 1-p
  mnormp <- cbind(p,1-p)
  M <- 2
  #  convert to surprisal values
  mnorms <- -log2(mnormp)
  #  plot the surprisal values
  matplot(x, mnorms, type="l", lty=c(1,1), col=c(1,1), 
          ylab="Surprisal (2-bits)")
  # add some log-normal error
  mnormdata <- exp(log(mnorms) + rnorm(11)*0.1)
  #  set up a b-spline basis object
  nbasis <- 7
  sbasis <- create.bspline.basis(c(-2,2),nbasis)
  #  define an initial coefficient matrix
  cmat <- matrix(0,7,1)
  #  set up a fd object for suprisal smoothing
  Sfd <- fd(cmat, sbasis)
  Zmat <- matrix(c(1,-1),2,1)
  #  smooth the noisy data
  result <- smooth.surp(x, mnormdata, cmat, Sfd, Zmat)
  #  plot the data and the fits of the two surprisal curves
  xfine <- seq(-2,2,len=51)
  sfine <- eval.surp(xfine, result$Sfd, Zmat)
  matplot(xfine, sfine, type="l", lty=c(1,1), col=c(1,1))
  points(x, mnormdata[,1])
  points(x, mnormdata[,2])
  #  convert the surprisal fit values to probabilities
  pfine <- 2^(-sfine)
  #  check that they sum to one
  apply(pfine,1,sum)
  par(oldpar)
}

