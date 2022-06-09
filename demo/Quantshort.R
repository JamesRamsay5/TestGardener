
#  -------------  preliminaries  -----------------

library(TestGardener)
#  set working directory to TestGardener

# ----------- read in data and key ----------- 

titlestr  <- "SweSAT-Q: 24 math analysis items"

U         <- scan("data/Ushort.txt", "o") # used from TestGardener/R
N         <- length(U) # Number of examinees
Umat      <- as.integer(unlist(stringr::str_split(U,"")))
n         <- length(Umat)/N # Number of items
U         <- matrix(Umat,N,n,byrow=TRUE)

Quantshort_U <- U

save(Quantshort_U, file="data/Quantshort_U.rda")  #  17 May 2022
load(Quantshort_U, file="data/Quantshort_U.rda")

knostring <- scan("data/keyshort.txt", "o")
key       <- as.integer(unlist(stringr::str_split(knostring,"")))

Quantshort_key <- key

save(Quantshort_key, file="data/Quantshort_key.rda")  #  17 May 2022
load(Quantshort_key, file="data/Quantshort_key.rda")

# Define the coding of invalid responses.
# Any value in U greater than number of legitimate options is treated as missing.

noption <- rep(0,n)
for (i in 1:n)
{
  noption[i]  <- length(unique(U[,i]))
}

grbg <- noption  # assumed that the last option is always a garbage option

# summary count for each option

for (i in 1:n)
{
  print(paste("Item: ",i,sep=""))
  for (m in 1:noption[i])
  {
    print(paste("Option: ",m," N= ",sum(U[,i]==m),sep=""))
  }
}

# --------- Define the option score values for each item ---------

ScoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  ScoreList[[item]] <- scorei
}

optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

# ----------------  Initialization Steps  ------------------------

Quantshort_dataList <- make.dataList(U, key, optList, grbg, titlestr=titlestr)

# Quantshort_dataList <- make.dataList(U, key, optList, grbg, scrrng, titlestr, 
#                                 nbin, Wnbasis, jitterwrd=FALSE, 
#                                 linearwrd=TRUE)

#  save this list object in the data folder

save(Quantshort_dataList, file="data/Quantshort_dataList.rda")  #  17 May 2022
load(file="data/Quantshort_dataList.rda")

#  plot the sum scores as a histogram 

hist(Quantshort_dataList$scrvec, Quantshort_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

#  compute the initial option surprisal curves using the 
#  percentage ranks as initial estimates of theta

theta     <- Quantshort_dataList$percntrnk
thetaQnt  <- Quantshort_dataList$thetaQnt

WfdResult <- Wbinsmth(theta, Quantshort_dataList)

#  Plot the initial option proability and surprisal curves

WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)

indfine   <- seq(0,100,len=101)
plotindex <- 1:n
plotrange <- c(0,100)

Wbinsmth.plot(indfine, plotindex, plotrange, binctr, Qvec, Quantshort_dataList, WfdList, 
              Wrng=c(0,3))
ICC.plot(indfine, WfdList, Quantshort_dataList, Qvec, binctr,  Wrng=c(0,3))

# ---------------  Optimal scoring: cycle of smoothing/theta estimation  ------------

#  Set number of cycles and the cell array to containing the parameter

ncycle=10

#  ----------------------------------------------------------------------------
#                      Proceed through the cycles
#  ----------------------------------------------------------------------------

AnalyzeResult <- Analyze(theta, thetaQnt, Quantshort_dataList, ncycle, 
                         itdisp=FALSE, verbose=FALSE) 

parList  <- AnalyzeResult$parList
meanHvec <- AnalyzeResult$meanHsave
  
#  ----------------------------------------------------------------------------
#              Plot the average H value, meanHsave, over cycles
#  ----------------------------------------------------------------------------

cycleno <- 1:ncycle
plot(cycleno,meanHvec[cycleno], type="b", lwd=2, xlab="Cycle Number")

#  select cycle for plotting

icycle <- 10

Quantshort_parList  <- parList[[icycle]]

#  save this list object in the data folder

save(Quantshort_parList, file="data/Quantshort_parList.rda") # May 17, 2022
load(file="data/Quantshort_parList.rda")

WfdList    <- Quantshort_parList$WfdList
Qvec       <- Quantshort_parList$Qvec
binctr     <- Quantshort_parList$binctr
theta      <- Quantshort_parList$theta

#  ----------------------------------------------------------------------------
#                   Compute the arc length or information measure 
#  ----------------------------------------------------------------------------

# A variety of useful objects related to the test info manifold
# are computed and returned in a struct object infoStr

Quantshort_infoList <- theta2arclen(theta, Qvec, WfdList, binctr)

save(Quantshort_infoList, file="data/Quantshort_infoList.rda")  # 20 May 2022
load(file="data/QuantShort_infoList.rda")

#  The length of the test manifold
arclength     <- Quantshort_infoList$arclength
#  The log derivative fd object for calculating arc length values from 
#  thetavalues
Wfd_theta     <- Quantshort_infoList$Wfd_theta
#  indefinite integral of arc length values corresponding to 
#  equally spaced theta values
arclengthvec  <- Quantshort_infoList$arclengthvec
#  The N arc length values corresponding to the N estimated score indes
#  theta values
theta_al      <- Quantshort_infoList$theta_al
#  The arc length values for the five marker percentages
Qvec_al       <- Quantshort_infoList$Qvec_al
#  The arc length values for the five marker percentages
binctr_al       <- Quantshort_infoList$binctr_al
#  The log derivative fd object for calculating theta values from 
#  arclength values
Wfd_info      <- Quantshort_infoList$Wfd_info
#  101 score index values corresponding to 101 equally spaced arc length
#  values
thetavec      <- Quantshort_infoList$thetavec
#  The dimension of the overspace within which the test info manifold
#  is found.
Wdim          <- Quantshort_infoList$Wdim

#  ----------------------------------------------------------------------------
#                   Plot surprisal curves for each test question
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

#  plot both the probability and surprisal curves along with data points

indfine   <- seq(0,100,len=101)
plotindex <- 1:n
plotrange <- c(0,100)

#  over score index theta
ICC.plot(indfine, WfdList, Quantshort_dataList, Qvec, binctr,  Wrng=c(0,5))

#  over arclength or information
ICC.plot(arclengthvec, WfdList, Quantshort_dataList, Qvec_al, binctr_al,  Wrng=c(0,2.5))

#  ----------------------------------------------------------------------------
#                         Plot density of theta
#  ----------------------------------------------------------------------------

ttllab     <- paste(titlestr,": percent rank", sep="")
scrrng     <- c(0,100)
indden10   <- scoreDensity(theta, scrrng, ttlstr=ttllab)
print(indden10)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu <- testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- scoreDensity(mu, Quantshort_dataList$scrrng, ttlstr=ttllab) 
print(muden)

#  compute expected score for each value in the fine mesh of theta values

indfine <- seq(0,100,len=101)
mufine  <- testscore(indfine, WfdList, optList)
mu.plot(mufine, Quantshort_dataList$scrrng, titlestr)

#  ----------------------------------------------------------------------------
#         Compute arc length over a fine mesh of theta values and plot
#  ----------------------------------------------------------------------------

#  print length of the test effort curve

print(paste("Arc length =", round(arclength,2)))

#  plot arc length over fine mesh

ArcLength.plot(arclength, arclengthvec, titlestr)

#  ----------------------------------------------------------------------------
#         plot the distribution of score index values theta 
#  ----------------------------------------------------------------------------

#  Use a histogram to see fine detail

hist(theta,51)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
             titlestr="SweSAT 13B Theta Density",  
             scrnbasis=11, nfine=101)

#  ----------------------------------------------------------------------------
#         plot the distribution of arclength or information
#  ----------------------------------------------------------------------------

#  Use a histogram to see fine detail

hist(theta_al,51)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
             titlestr="SweSAT 13B Info Density",  
             scrnbasis=15, nfine=101)
#  ----------------------------------------------------------------------------
#  Display test effort curve projected into its first two principal components
#  ----------------------------------------------------------------------------

# nharm=2
Result <- Wpca.plot(arclength, WfdList, Quantshort_dataList$Wdim, titlestr=titlestr)

# nharm=3
Result <- Wpca.plot(arclength, WfdList, Quantshort_dataList$Wdim, 3, dodge = 1.005, 
                    titlestr=titlestr)

#  ----------------------------------------------------------------------------
#                          Display sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

Sensitivity.plot(arclengthvec, WfdList, Qvec_al, Quantshort_dataList)
  
#  ----------------------------------------------------------------------------
#                          Display power curves
#  ----------------------------------------------------------------------------

Power.plot(arclengthvec, WfdList, Qvec_al, Quantshort_dataList, height=0.25)

#  ----------------------------------------------------------------------------
#                          Display entropy curves
#  ----------------------------------------------------------------------------

Entropy.plot(arclengthvec, WfdList, Qvec_al, Quantshort_dataList, height=1)

#  ----------------------------------------------------------------------------
#                   Display H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

index <- 1:5

Hfuns.plot(theta, WfdList, Quantshort_dataList$U, plotindex=1:5)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

simList <- dataSimulation(Quantshort_dataList, Quantshort_parList, nsample=500)
dataSimulation.plot(simList, Qvec)
