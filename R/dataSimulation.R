dataSimulation <- function(dataList, parmList, nsample=1000) {
  
  #  Simulate data, analyze each simulated sample, and print results
  
  #  Last modified on 16 November 2023
  
  #  --------------------------------------------------------------------------
  #                   Set up access to population info:
  #  --------------------------------------------------------------------------
  
  #  info from dataList
  
  ScoreList <- dataList$ScoreList
  scrrng    <- dataList$scrrng
  Wdim      <- dataList$Wdim
  
  #  info from parmList
  
  SfdList   <- parmList$SfdList
  logdensfd <- parmList$logdensfd
  Qvec      <- parmList$Qvec
  binctr    <- parmList$binctr
  n         <- length(SfdList)
  
  #  --------------------------------------------------------------------------
  #                       Generate simulated data
  #  --------------------------------------------------------------------------
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  #  define population index values 
  
  nindex <- nfine
  
  #  compute population values for index, mu and arc length
  #  corresponding to nfine equally spaced values of index.
  
  index.pop=seq(0,100,len=101)
  mu.pop       <- mu(index.pop, SfdList, dataList$optList)
  result       <- index2info(index.pop, Qvec, SfdList, binctr)
  scopevec.pop <- result$scopevec 
  infoSurp     <- result$infoSurp 
  infoSurp1001 <- result$infoSurp1001 
  infoSurpfine <- seq(0,infoSurp,len=101)
  
  #  arrays to store simulated data
  
  chcematsave <- array(0,c(nindex,n,nsample))
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) print(paste("Sample ",isample))
    chcematsave[,,isample] <- chce_simulate(index.pop, SfdList)
  }
  
  #  --------------------------------------------------------------------------
  #                    Analyze each sample
  #  --------------------------------------------------------------------------
  
  #  Note that total arc length is a function of the surprisal functions,
  #  and therefore is fixed for these simulations.
  
  #  set up matrices to save results
  
  sumscrsave <- matrix(0, nindex, nsample)
  indexsave  <- matrix(0, nindex, nsample)
  musave     <- matrix(0, nindex, nsample)
  infosave   <- matrix(0, nindex, nsample)
  
  Fvalsave   <- matrix(0, nindex, nsample)
  DFvalsave  <- matrix(0, nindex, nsample)
  D2Fvalsave <- matrix(0, nindex, nsample)
  
  Infoopt <- FALSE  #  Use info matrix instead of Hessian
  
  #  analyze the samples
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) print(paste("Sample ", isample))
    # print(isample)
    chcemati <- chcematsave[,,isample]
    #  sum scores for this sample
    # scrveci <- sumscorefn(chcemati, ScoreList)
    scrveci <- matrix(0,nindex,1)
    for (j in 1:nindex) {
      for (i in 1:n) {
        scrveci[j] = scrveci[j] + ScoreList[[i]][chcemati[j,i]]
      }
    }
    #  optimal scores for multi-option analysis
    result     <- index_fun(index.pop, SfdList, chcemati)             
    indexveci  <- result$index_out 
    Fvalveci   <- result$Fval 
    DFvalveci  <- result$DFval 
    D2Fvalveci <- result$D2Fval
    muveci     <- mu(indexveci, SfdList, dataList$scoreList)
    infoveci   <- pracma::interp1(as.numeric(indfine), 
                                  as.numeric(infoSurpfine), 
                                  as.numeric(indexveci))
    sumscrsave[,isample] <-    scrveci
    indexsave[,isample]  <-  indexveci     
    musave[,isample]     <-     muveci
    infosave[,isample]   <-   infoveci
    Fvalsave[,isample]   <-   Fvalveci
    DFvalsave[,isample]  <-  DFvalveci
    D2Fvalsave[,isample] <- D2Fvalveci
  }
  
  #  --------------------------------------------------------------------------
  #             Define list  object simList to save the results
  #  --------------------------------------------------------------------------
  
  simList <- list(
    sumscr     = sumscrsave,
    index      = indexsave,
    mu         = musave,
    info       = infosave,
    index.pop  = index.pop,
    mu.pop     = mu.pop,
    n          = n,
    nindex     = nindex,
    indfine    = indfine, 
    Qvec       = Qvec
  )
  
  return(simList)
  
}
