#  -------------  preliminaries  -----------------

library(TestGardener)
#  set working directory to TestGardener

# ----------- read in data ----------- 

titlestr  <- "Symptom Distress Scale"

U         <- scan("SDS_U.txt", "o") # used from TestGardener/R
U         <- matrix(U,473,13,byrow=TRUE)

# U         <- U[,2]
# N         <- length(U) # Number of examinees
# Umat      <- as.integer(unlist(stringr::str_split(U,"")))
# n         <- length(Umat)/N # Number of items
# U         <- matrix(Umat,N,n,byrow=TRUE)

key     <- NULL

noption <- matrix(5,n,1)

# Change the coding of invalid responses.
for (i in 1:n)
{
  if (any(U[,i] > noption[i]))
  {
    noption[i]  <- noption[i] + 1 # Add one option for invalid responses
    U[U[,i] >= noption[i],i] <- noption[i]
  }
}

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
  scorei <- c(0:4,0)
  ScoreList[[item]] <- scorei
}

optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

# ----------------  Initialization Steps  ------------------------

SDS_dataList <- make.dataList(U, key, optList, titlestr=titlestr)

#  save this list object in the data folder

save(SDS_dataList, file="data/SDS_dataList.rda")
load(file="data/SDS_dataList.rda")

#  --------- Set initial values that are required in the later analysis --------- 

#  compute the initial option surprisal curves using the 
#  percentage ranks as initial estimates of theta

theta     <- SDS_dataList$percntrnk
thetaQnt  <- SDS_dataList$thetaQnt

WfdResult <- Wbinsmth(theta, SDS_dataList)

#  Plot the initial option proability and surprisal curves

WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)

indfine   <- seq(0,100,len=101)
plotindex <- 1:n
plotrange <- c(0,100)

Wbinsmth.plot(indfine, plotindex, plotrange, binctr, Qvec, Quant_dataList, WfdList, 
              Wrng=c(0,2.5))
ICC.plot(indfine, WfdList, SDS_dataList, Qvec, binctr,  Wrng=c(0,2.5))

#  re-display the curves for the first scale item

print(plot_list[[1]])

# ---------------  Optimal scoring: cycle of smoothing/theta estimation  ------------

#  Set number of cycles and the cell array to containing the parameter

ncycle=20

#  ----------------------------------------------------------------------------
#                      Proceed through the cycles
#  ----------------------------------------------------------------------------

AnalyzeResult <- Analyze(theta, thetaQnt, SDS_dataList, ncycle, itdisp=FALSE) 

parList   <- AnalyzeResult$parList
meanHsave <- AnalyzeResult$meanHsave

#  ----------------------------------------------------------------------------
#              Plot meanHsave and choose cycle for plotting
#  ----------------------------------------------------------------------------

cycleno <- 1:ncycle
par(mfrow=c(1,1))
plot(cycleno,meanHsave, type="b", lwd=2, xlab="Cycle Number")

#  select cycle for plotting

icycle <- 20

SDS_parList  <- parList[[icycle]]

#  save this list object in the data folder

save(SDS_parList, file="data/SDS_parList.rda")  #  17 May 2022
load(file="data/SDS_parList.rda")  

WfdList    <- SDS_parList$WfdList
Qvec       <- SDS_parList$Qvec
binctr     <- SDS_parList$binctr
theta      <- SDS_parList$theta

#  ----------------------------------------------------------------------------
#                   Compute the arc length or information measure 
#  ----------------------------------------------------------------------------

# A variety of useful objects related to the test info manifold
# are computed and returned in a struct object infoStr

SDS_infoList <- theta2arclen(theta, Qvec, WfdList)  # saved 6 May 2022

save(SDS_infoList, file="data/SDS_infoList.rda")  # 17 May 2022
load(file="SDS_infoList.rda")

#  The length of the test manifold
arclength     <- SDS_infoList$arclength
#  The log derivative fd object for calculating arc length values from 
#  thetavalues
Wfd_theta     <- SDS_infoList$Wfd_theta
#  indefinite integral of arc length values corresponding to 
#  equally spaced theta values
arclengthvec  <- SDS_infoList$arclengthvec
#  The N arc length values corresponding to the N estimated score indes
#  theta values
theta_al      <- SDS_infoList$theta_al
#  The arc length values for the five marker percentages
Qvec_al       <- SDS_infoList$Qvec_al
#  The arc length values for the bin centers
binctr_al     <- SDS_infoList$binctr_al
#  The log derivative fd object for calculating theta values from 
#  arclength values
Wfd_info      <- SDS_infoList$Wfd_info
#  101 score index values corresponding to 101 equally spaced arc length
#  values
thetavec      <- SDS_infoList$thetavec
#  The dimension of the overspace within which the test info manifold
#  is found.
Wdim          <- SDS_infoList$Wdim

#  ----------------------------------------------------------------------------
#                   Plot surprisal curves for each test question
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

#  over score index theta

indfine   <- seq(0,100,len=101)
ICC.plot(indfine, WfdList, SDS_dataList, Qvec, binctr,  Wrng=c(0,2.5))

#  over arclength or information

ICC.plot(arclengthvec, WfdList, SDS_dataList, Qvec_al, binctr_al,  
         Wrng=c(0,2.5))

#  ----------------------------------------------------------------------------
#                         Plot density of theta
#  ----------------------------------------------------------------------------

ttllab     <- paste(titlestr,": percent rank", sep="")
edges      <- c(0,100)
theta_in   <- theta[0 < theta & theta < 100]
indden10   <- scoreDensity(theta_in, edges, 15, ttlstr=ttllab)
print(indden10)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu <- testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
scrrng <- c(0,37)
muden  <- scoreDensity(mu, scrrng, ttlstr=ttllab) 
print(muden)

#  compute expected score for each value in the fine mesh of theta values

indfine <- seq(0,100,len=101)
mufine  <- testscore(indfine, WfdList, optList)
mu.plot(mufine, SDS_dataList$scrrng, ttllab)

#  ----------------------------------------------------------------------------
#         Compute arc length over a fine mesh of theta values and plot
#  ----------------------------------------------------------------------------

#  print length of the test information curve

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
             titlestr="SDS Theta Density",  
             scrnbasis=11, nfine=101)

#  ----------------------------------------------------------------------------
#         plot the distribution of arclength or information
#  ----------------------------------------------------------------------------

#  Use a histogram to see fine detail

hist(theta_al,51)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
             titlestr="SDS Info Density",  
             scrnbasis=15, nfine=101)

#  ----------------------------------------------------------------------------
#  Display test effort curve projected into its first two principal components
#  ----------------------------------------------------------------------------

# nharm=2
Result <- Wpca.plot(arclength, WfdList, SDS_dataList$Wdim, titlestr=titlestr)

# nharm=3
Result <- Wpca.plot(arclength, WfdList, SDS_dataList$Wdim, 3, dodge = 1.005, 
                    titlestr=titlestr)

#  ----------------------------------------------------------------------------
#                          Display sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

#  over score index theta

Sensitivity.plot(indfine, WfdList, Qvec, SDS_dataList, 
                 titlestr=titlestr, plotindex=1:n)

#  over arclength or information

Sensitivity.plot(arclengthvec, WfdList, Qvec_al, SDS_dataList, 
                 titlestr=titlestr, plotindex=1:n, plotrange=c(0,arclength))

#  ----------------------------------------------------------------------------
#                          Display power curves
#  ----------------------------------------------------------------------------

#  over score index theta

Power.plot(indfine, WfdList, Qvec, SDS_dataList, plotindex=1:n, height=0.3)

#  over arclength or information

Power.plot(arclengthvec, WfdList, Qvec_al, SDS_dataList, plotindex=1:n, 
           plotrange=c(0,arclength), height=0.3)

#  ----------------------------------------------------------------------------
#                          Display entropy curves
#  ----------------------------------------------------------------------------

#  over score index theta

Entropy.plot(indfine, WfdList, Qvec, SDS_dataList, height=1)

#  over arclength or information

Entropy.plot(arclengthvec, WfdList, Qvec_al, SDS_dataList, plotindex=1:n, 
             plotrange=c(0,arclength), height=1)

#  ----------------------------------------------------------------------------
#                   Display H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

Hfuns.plot(theta, WfdList, SDS_dataList$U, plotindex=1:5)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

simList <- dataSimulation(SDS_dataList, SDS_parList, nsample=500)
dataSimulation.plot(simList, Qvec)

