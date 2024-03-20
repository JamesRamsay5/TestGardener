surprisal.chart <- function(B, Zmat, C = NULL) {
  #  This is a retraction function for mapping real matrix y of M-1
  #  columns and N rows on to the probability (P) and surprisal (S) 
  #  manifolds within M-space.
  #  x = y*Z' is a real vector of length M where Z is an M by M-1 matrix
  #  with orthonormal columns that sum to zero.  
  #  Function zerobasis(M) generates Z using fourier basis values.
  #  Vectors in N by M matrix X are transformed into
  #  a multinomial vector containing probability M values that sum to one, 
  #  and  also into surprisal, S = -log_C P.
  #  C is the base for the surprisal transform, defaults to M.
  #  S is the surprisal vector.
  #  P is the probability vector, that is, a multinomial object.
  #  DS is the matrix of partial derivatives of surprisal vector S elements
  #     with respect to elements of the real vector x.
  
  # Last modified 2023-02-24 by Juan Li
  # Based on the Matlab code (2023-01-10)
  
  N <- nrow(B)
  M <- nrow(Zmat)
  if (ncol(Zmat) != M-1)
  {
    stop(paste("Zmat does not have M-1 columns: M = ", M, 
               " ncol(Zmat) = ", ncol(Zmat), sep = ""))
  }
  
  if (ncol(B) != M-1)
  {
    stop(paste("B does not have M-1 columns: M = ", M, 
               " ncol(B) = ", ncol(B), sep = ""))
  }
  
  if (is.null(C)) C <- M
  logC <- log(C)
  
  # convert y into x
  X <- B %*% t(Zmat)

  # compute probability and surprisal vectors

  E <- C^(-X)
  sumE <- rowSums(E)
  matsumE <- matrix(sumE, nrow=length(sumE), ncol=M, byrow=FALSE)
  P <- E/matsumE
  S <- -log(P)/logC

  #  compute, if needed, the partial derivatives of S with respect to x
  DS <- array(rep(0, N*M*M), dim = c(N, M, M))
  
  for (i in 1:M)
  {
    for (j in 1:M)
    {
      if (i == j)
      {
        for (k in 1:N)
        {
          DS[k, j, i] <- 1 - P[k, j]
        }
      } else
      {
        for (k in 1:N)
        {
          DS[k, j, i] <- -P[k, j]
        }
      }
    }
  }
  
  return(list(S = S, P = P, DS = DS, X = X))
}