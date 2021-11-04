growsurf = function(tvec, a, b, c) {
  #  compute residual covariance surface resembling that for
  #  growth residuals at time points in TVEC.
  #  Last modified 29 September 2011

  n = length(tvec)
  sigmat = matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      spt = (tvec[i] + tvec[j])/2
      smt = (abs(tvec[i]-tvec[j]))
      sigmat[i,j] = (1+2*exp(-spt/a) +
                       0.5*exp(-(spt-9)^2/4))*cos(2*pi*smt/b)*exp(-(smt/c))/4
    }
  }
  return(sigmat)
}
