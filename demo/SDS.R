#  -------------  preliminaries  -----------------

library(TestGardener)
#  set working directory to TestGardener

# ----------- read in data ----------- 

titlestr  <- "Symptom Distress Scale"
U         <- scan("../data/SDS.txt", "o") # used from TestGardener/R
# U         <- scan("inst/SDS.txt", "o")
U         <- matrix(U,473,2,byrow=TRUE)
U         <- U[,2]
N         <- length(U) # Number of examinees
Umat      <- as.integer(unlist(stringr::str_split(U,"")))
n         <- length(Umat)/N # Number of items
U         <- matrix(Umat,N,n,byrow=TRUE)

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

optList <- list() # option scores
for (item in 1:n){
  scorei <- c(0:4,0)
  optList[[item]] <- scorei
}

optList <- list(itemLab=NULL, optLab=NULL, optScr=optList)

SDS_dataList <- make.dataList(U, key, optList, grbg, scrrng = c(0,37))

#  save this list object in the data folder

save(SDS_dataList, file="data/SDS_dataList.rda")

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

plot_list <- Wbinsmth.plot(binctr, Qvec, WfdList, SDS_dataList, Wrng=c(0,3))

#  re-display the curves for the first scale item

print(plot_list[[1]])

# ---------------  Optimal scoring: cycle of smoothing/theta estimation  ------------

#  Set number of cycles and the cell array to containing the parameter

ncycle=10

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

icycle <- 10

SDS_parList  <- parList[[icycle]]

WfdList    <- SDS_parList$WfdList
Qvec       <- SDS_parList$Qvec
binctr     <- SDS_parList$binctr
theta      <- SDS_parList$theta
arclength  <- SDS_parList$arclength
alfine     <- SDS_parList$alfine

#  save this list object in the data folder

save(SDS_parList, file="../data/SDS_parList.rda")

#  ----------------------------------------------------------------------------
#                   Plot surprisal curves for each test question
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

Wmax <- 2.5
Wbinsmth.plot(binctr, Qvec, WfdList, SDS_dataList, Wrng=c(0,Wmax))

#  ----------------------------------------------------------------------------
#                         Plot density of theta
#  ----------------------------------------------------------------------------

ttllab     <- paste(titlestr,": percent rank", sep="")
edges      <- c(0,100)
theta_in   <- theta[0 < theta & theta < 100]
indden10   <- scoreDensity(theta_in, edges, 15, ttlstr=ttllab)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu <- testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
scrrng <- c(0,37)
muden  <- scoreDensity(mu, scrrng, ttlstr=ttllab) 

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

ArcLength.plot(arclength, alfine, titlestr)

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

Sensitivity.plot(WfdList, Qvec, SDS_dataList, titlestr=titlestr, plotindex=1:n)

#  ----------------------------------------------------------------------------
#                          Display power curves
#  ----------------------------------------------------------------------------

Power.plot(WfdList, Qvec, SDS_dataList, plotindex=1:n, height=0.3)

#  ----------------------------------------------------------------------------
#                   Display H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

Hfuns.plot(theta, WfdList, SDS_dataList$U, plotindex=1:5)

