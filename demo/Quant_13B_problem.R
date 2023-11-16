#  Data consisting only of the 24 math problem  solving items:  c(1:12,41:52)
#  Used to generate the vignette version May 4, 2023

#  -------------  preliminaries  -----------------

titlestr  <- "SweSAT-Q_13B_problem"

#  --------------  starting from scratch  -------------------------------------
#. both the choice indices in matrix U and the right answer key
#. are stored as a string of  integers with no spaces
#. The stringr package is used to convert the strings to 
#. integer vectors

chcevec <- as.numeric(scan("data/chcemat.txt","o"))
chcemat <- matrix(chcevec,1000,24,byrow=TRUE)
N <- nrow(chcemat)
n <- ncol(chcemat)
knostring <- scan("data/key.txt", "o") 
key       <- as.integer(unlist(stringr::str_split(knostring,"")))
noption <- rep(0,n)
for (i in 1:n) noption[i]  <- 4
scoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  scoreList[[item]] <- scorei
}

TGresult <- TGanalysis(chcemat, scoreList, noption, 
                       NumBasis=4, ncycle=10, verbose=TRUE)

dataList <- TGresult$dataList
parmList <- TGresult$parmList
infoList <- TGresult$infoList

SfdList     <- parmList$SfdList
Qvec        <- parmList$Qvec
binctr      <- parmList$binctr
indxScr     <- parmList$indxScr
Hval        <- parmList$Hval
DHval       <- parmList$DHval
D2Hval      <- parmList$D2Hval

infoSurp    <- infoList$infoSurp
infoSurpvec <- infoList$infoSurpvec
scopevec    <- infoList$scopevec
Qinfovec    <- infoList$Qinfovec
bininfoctr  <- infoList$bininfoctr

#  ----------------------------------------------------------------------------
#        Plot surprisal curves for each test question using ICC_plot
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

indfine   <- seq(0,100,len=101)
plotindex <- 1:n
plotrange <- c(0,100)

#  over score index indxScr

plotindex=6
ICC_plot(indfine, SfdList, dataList, Qvec, binctr, 
         data_point=TRUE, plotType=c("P","S"), Srng=c(0,4), 
         DSrng=c(-0.2, 0.2), plotindex=plotindex)

#  over infoSurp or information

ICC_plot(infoSurpvec, SfdList, dataList, Qinfovec, bininfoctr,  
         data_point=TRUE, plotType=c("P","S"), Srng=c(0,4))

#  ----------------------------------------------------------------------------
#               Plot density of indxScr  using scoreDensity
#  ----------------------------------------------------------------------------

ttllab   <- paste(titlestr,": percent rank", sep="")
scrrng   <- c(0,100)
indden10 <- scoreDensity(indxScr, scrrng, ttlstr=ttllab)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu     <- testscore(indxScr, SfdList, dataList$optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- scoreDensity(mu, dataList$scrrng, ttlstr=ttllab) 

#  compute expected score for each value in the fine mesh of indxScr values

indfine <- seq(0,100,len=101)
mufine  <- testscore(indfine, SfdList, dataList$optList)
mu_plot(mufine, dataList$scrrng, titlestr)

#  ----------------------------------------------------------------------------
#         Plot arc length over a fine mesh of indxScr values and plot
#  ----------------------------------------------------------------------------

#  print length of the test effort curve

print(paste("Arc length =", round(infoSurp,2)))

#  plot arc length over fine mesh

ArcLength_plot(infoSurp, infoSurpvec, titlestr)

#  ----------------------------------------------------------------------------
#         plot the distribution of score index values indxScr 
#  ----------------------------------------------------------------------------

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

par(mfrow=c(2,1))
density_plot(indxScr, c(0,100), Qvec, xlabstr="Score index", 
             titlestr="SweSAT 13B Theta Density",  
             scrnbasis=11, nfine=101)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(scopevec, c(0,infoSurp), Qinfovec, xlabstr="Arclength", 
             titlestr="SweSAT 13B Info Density",  
             scrnbasis=15, nfine=101)

#  ----------------------------------------------------------------------------
#  print test info curve from its first two or three principal components
#  ----------------------------------------------------------------------------

# nharm=2
Results <- Spca_plot(infoSurp, SfdList, dataList$Sdim, 
                     titlestr=titlestr)
print(paste("Variance percent = ",100*round(Results[[3]],3)))

# nharm=3
Results <-Spca_plot(infoSurp, SfdList, dataList$Sdim,   
                    3, dodge = 1.005, titlestr=titlestr)
print(paste("Variance percent = ",100*round(Results[[3]],3)))

#  ----------------------------------------------------------------------------
#                          print sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

Sensitivity_plot(indfine, SfdList, Qvec, dataList)
  
#  ----------------------------------------------------------------------------
#                          print power curves
#  ----------------------------------------------------------------------------

Power_plot(indfine, SfdList, Qvec, dataList, height=0.25)

#  ----------------------------------------------------------------------------
#                          print entropy curves
#  ----------------------------------------------------------------------------

Entropy_plot(indfine, SfdList, Qvec, dataList, height=1)

#  ----------------------------------------------------------------------------
#                   printlay H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

Ffuns_plot(indfine, theta, SfdList, dataList, plotindex=1:5)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

# simList <- dataSimulation(dataList, parmList, nsample=50)
# 
# scorePerformance(simList, Qvec)
# 
# ##  Item analysis plots over peak intervals
# 
# peakbdry <- C(0,21,47,67,100)
# 
# Srng <- C(0,5)

#  ----------------------------------------------------------------------------
#                Compute item arc lengths over [0,100]  
#  ----------------------------------------------------------------------------

itemArclength_0_100  <- matrix(0,80,1)
indfine <- seq(0,100,len=1000)
indxScri <- NULL
for (i in 1:n) {
  print(i)
  Resulti   <- indxScr2arclen(indxScr, Qvec, SfdList, binctr, i)
  indxScri  <- Resulti$scopevec
  infoSurpi <- Resulti$infoSurp
  itemArclength_0_100[i] <- infoSurpi
}

#  print all of the arc lengths

print("infoSurp_0_100:")
print(itemArclength_0_100)

#  histogram of the 80 arc lengths 

hist(itemArclength_0_100)

#  sort 80 infoSurps in descending order

Result <- sort.int(itemArclength_0_100, decreasing=TRUE, index.return=TRUE)
alsort <- Result$x 
isort  <- Result$ix

#  printlay the four items with the longest arc lengths:
#  28, 63, 24, 26

#  plot the surprisal curves in descending order of arc length

Srng <- c(0,5)

ICC_plot(SfdList, dataList, Qvec, binctr, plotType = c("P", "S"), 
         plotindex=isort)

#  plot the surprisal curves in ascending order of arc length

Result <- sort.int(itemArclength_0_100)
alsort <- Result$alsort
isort  <- Resuslt$isort
ICC_plot(SfdList, dataList, Qvec, binctr, plotType = c("P", "S"), 
         plotindex=isort)

#  ----------------------------------------------------------------------------
#               Compute item arc lengths over the peak intervals 
#  ----------------------------------------------------------------------------

itemArclength_80_100 <- matrix(0,80,1)

indfine <- seq(80,100,len=101)

pltrng <- c(80,100)
indxScri <- NULL
for (i in 1:80) {
  Result     <- indxScr2arclen(indxScr, Qvec, SfdList, binctr, i, pltrng, TRUE)
  indxScri     <- Result$scopevec
  infoSurpi <- Result$infoSurp
  itemArclength_80_100[i] <- infoSurpi
}

# print('infoSurp_80_100:')
# print([c((1:80),itemArclength_80_100])

hist(itemArclength_80_100)

Result <- sort.int(itemArclength_80_100, decreasing=TRUE, index.return=TRUE)
alsort <- Result$x 
isort  <- Result$ix

ICC_plot(indfine, SfdList, dataList, Qvec, binctr, 
         plotType = c("P", "S"), plotindex=isort, plotrange=c(80,100))


