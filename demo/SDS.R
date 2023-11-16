#  -------------  preliminaries  -----------------

library(TestGardener)
#  set working directory to TestGardener

# ----------- read in data ----------- 

titlestr  <- "Symptom Distress Scale"

chcevec <- as.integer(scan("data/SDS_chcemat.txt", "o")) # used from TestGardener/R
chcemat <- matrix(chcevec,473,13,byrow=TRUE)

N <- nrow(chcemat)
n <- ncol(chcemat)

key     <- NULL

noption <- rep(5,n)

# Change the coding of invalid responses.
for (i in 1:n) {
  if (any(chcemat[,i] > noption[i]))
  {
    noption[i]  <- noption[i] + 1 # Add one option for invalid responses
    chcemat[chcemat[,i] >= noption[i],i] <- noption[i]
  }
}

# summary count for each option
for (i in 1:n) {
  print(paste("Item: ",i,sep=""))
  for (m in 1:noption[i]) {
    print(paste("Option: ",m," N= ",sum(chcemat[,i]==m),sep=""))
  }
}

# --------- Define the option score values for each item ---------

ScoreList <- vector("list", n) # option scores
for (item in 1:n){
  scorei <- c(0:4,0)
  ScoreList[[item]] <- scorei
}

TGresult <- TGanalysis(chcemat, scoreList, noption, 
                       NumBasis=4, ncycle=10, verbose=TRUE)

SDS_dataList <- TGresult$dataList
SDS_parmList <- TGresult$parmList
SDS_infoList <- TGresult$infoList

#  save this list object in the data folder 5 May 2023

save(SDS_dataList, 
     file="../TestGardener/data/SDS_dataList.rda") 
load(file="../TestGardener/data/SDS_dataList.rda") 

save(SDS_parmList, 
     file="../TestGardener/data/SDS_parmList.rda") 
load(file="../TestGardener/data/SDS_parmList.rda") 

save(SDS_infoList, 
     file="../TestGardener/data/SDS_infoList.rda")  
load(file="../TestGardener/data/SDS_infoList.rda")

SfdList    <- SDS_parmList$SfdList
Qvec       <- SDS_parmList$Qvec
binctr     <- SDS_parmList$binctr
index      <- SDS_parmList$index
Hval       <- SDS_parmList$Hval
DHval      <- SDS_parmList$DHval
D2Hval     <- SDS_parmList$D2Hval

infoSurp    <- SDS_infoList$infoSurp
infoSurpvec <- SDS_infoList$infoSurpvec
scopevec    <- SDS_infoList$scopevec
Qinfovec    <- SDS_infoList$Qinfovec
bininfoctr  <- SDS_infoList$bininfoctr

#  ----------------------------------------------------------------------------
#                   Plot surprisal curves for each test question
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

#  over score index index

indfine   <- seq(0,100,len=101)
ICC_plot(indfine, SfdList, SDS_dataList, Qvec, binctr,  Srng=c(0,2.5))

#  over infoSurp or information

ICC_plot(infoSurpvec, SfdList, SDS_dataList, Qinfovec, bininfoctr,  
         data_point=TRUE, plotType=c("P", "S"), Srng=c(0,2.5))

#  ----------------------------------------------------------------------------
#                         Plot density of index
#  ----------------------------------------------------------------------------

ttllab     <- paste(titlestr,": percent rank", sep="")
sccrrng    <- c(0,100)
index_int  <- index[0 < index & index < 100]
indden10   <- scoreDensity(index_int, scrrng=c(0,37), ttlstr=ttllab)
print(indden10)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

muvec <- mu(index, SfdList, SDS_dataList$scoreList)
ttllab <- paste(titlestr,": expected score", sep="")
scrrng <- c(0,37)
muden  <- scoreDensity(mu, scrrng, ttlstr=ttllab) 
print(muden)

#  compute expected score for each value in the fine mesh of index values

indfine <- seq(0,100,len=101)
mufine  <- testscore(indfine, SfdList, optList)
mu_plot(mufine, SDS_dataList$scrrng, ttllab)

#  ----------------------------------------------------------------------------
#         Compute arc length over a fine mesh of index values and plot
#  ----------------------------------------------------------------------------

#  print length of the test information curve

print(paste("Arc length =", round(infoSurp,2)))

#  plot arc length over fine mesh

Scope_plot(infoSurp, infoSurpvec, titlestr)

#  ----------------------------------------------------------------------------
#         plot the distribution of score index values index 
#  ----------------------------------------------------------------------------

#  Use a histogram to see fine detail

hist(index,51)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(index, c(0,100), Qvec, xlabstr="Score index", 
             titlestr="SDS Theta Density",  
             scrnbasis=11, nfine=101)

#  ----------------------------------------------------------------------------
#         plot the distribution of infoSurp or information
#  ----------------------------------------------------------------------------

#  Use a histogram to see fine detail

hist(infoSurpvec,51)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(infoSurpvec, c(0,infoSurp), Qinfovec, xlabstr="Arclength", 
             titlestr="SDS Info Density",  
             scrnbasis=15, nfine=101)

#  ----------------------------------------------------------------------------
#  Display test effort curve projected into its first two principal components
#  ----------------------------------------------------------------------------
#  Perform the singular value decomposition of the curve trajectory

ReturnPca3 <- Spca(SfdList, nharm=3, rotate=FALSE)

#  display the proportiions of variation and their sum

print(round(100*ReturnPca3$varpropvarmx,1))

#  96.5%  2.4%  0.6%

print(round(sum(100*ReturnPca3$varpropvarmx,1)))

#  100.0%

#  display the plot 

titlestr <- paste("SDS scale =", round(infoSurp,1),"(3-bits)")
harmvarmxfd <- ReturnPca3$harmvarmxfd


png("figs/ScaleInfo.png")
Spca_plot(harmvarmxfd, nharm=3, titlestr)
dev.off()

#  ----------------------------------------------------------------------------
#                          Display sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

#  over score index index

Sensitivity_plot(indfine, SfdList, Qvec, SDS_dataList, 
                 titlestr=titlestr, plotindex=1:n)

#  over infoSurp or information

Sensitivity_plot(infofine, SfdList, Qinfovec, SDS_dataList, 
                 titlestr=titlestr, plotindex=1:n, plotrange=c(0,infoSurp))

#  ----------------------------------------------------------------------------
#                          Display power curves
#  ----------------------------------------------------------------------------

#  over score index index

Power_plot(indfine, SfdList, Qvec, SDS_dataList, plotindex=1:n, height=0.3)

#  over infoSurp or information

Power_plot(infoSurpvec, SfdList, Qinfovec, SDS_dataList, plotindex=1:n, 
           plotrange=c(0,infoSurp), height=0.3)

#  ----------------------------------------------------------------------------
#                          Display entropy curves
#  ----------------------------------------------------------------------------

#  over score index index

Entropy_plot(indfine, SfdList, Qvec, SDS_dataList, height=1)

#  over infoSurp or information

Entropy_plot(scopevec, SfdList, Qinfovec, SDS_dataList, plotindex=1:n, 
             plotrange=c(0,infoSurp), height=1)

#  ----------------------------------------------------------------------------
#                   Display H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

Ffuns_plot(indfine, index, SfdList, SDS_dataList$chcemat, plotindex=11:20)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

simList <- dataSimulation(SDS_dataList, SDS_parmList, nsample=500)
dataSimulation_plot(simList, Qvec)

