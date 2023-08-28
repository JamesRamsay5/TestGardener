#  Data consisting only of the 24 math problem  solving items:  c(1:12,41:52)
#  Used to generate the vignette version May 4, 2023

#  -------------  preliminaries  -----------------

library(TestGardener)

#  set working directory to TestGardener_work/Quant

setwd("../Quant/Quant_13B_problem")

# ----------- read in data and key ----------- 

titlestr  <- "SweSAT-Q_13B_problem"

#  --------------  starting from scratch  -------------------------------------
#. both the choice indices in matrix U and the right answer key
#. are stored as a string of  integers with no spaces
#. The stringr package is used to convert the strings to 
#. integer vectors

U         <- scan("../Quant_13B_full.txt", "o") 
N         <- length(U) # Number of examinees
Umat      <- as.integer(unlist(stringr::str_split(U,"")))
n         <- length(Umat)/N # Number of items
U         <- matrix(Umat,N,n,byrow=TRUE)

#  Remove the examinees that took the test abroad

ind    = 266:N
U      = U[ind,]
N      = N - 265

#  change 8 and 9 into 5 and 6, respectively

for (j in 1:N) {
  for (i in c(23:28,63:68)) {
    if (U[j,i] == 8 || U[j,i] == 9) U[j,i] = 6
  }
  for (i in c(1:22,29:63,69:80)) {
    if (U[j,i] == 8 || U[j,i] == 9) U[j,i] = 5
  }
}

#   choose a random subsample of 1000 examinees and columns of
#   of math problem questions

index <- sample(1:N, 1000)
problemindex <- c(1:12,41:52)
U   <- U[index,problemindex]
N <- nrow(U)
n <- ncol(U)

#. the choice data are saved as a .rda file
#. named Quant_13B_problem_U.rda

save(U, file="../TestGardener/data/Quant_13B_problem_U.rda")  
load(file="../TestGardener/data/Quant_13B_problem_U.rda")

#. The choice index data are saved as a .txt file containing positive
#. integers separated by single spaces.  This file is used
#. in the vignette Quant_13B_problem

sink(file="Quant_13B_problem_U.txt")
for (j in 1:N) {
  for (i in 1:n) {
    cat(" ")
    cat(U[j,i])
  }
  cat("\n")
}
sink()

U <- as.matrix(read.table("Quant_13B_problem_U.txt"))

#. The choice data are input from the .txt file
#. named Quant_13B_problem_U.txt


#  ---------------------  input the key  -----------------------------------

knostring <- scan("Quant_13B_key.txt", "o") 
key       <- as.integer(unlist(stringr::str_split(knostring,"")))
key       <- key[problemindex]

save(key, file="../TestGardener/data/Quant_13B_problem_key.rda")  
load(file="../TestGardener/data/Quant_13B_problem_key.rda")

#. The key index data are saved as a .txt file containing positive
#. integers separated by single spaces.  This file is used
#. in the vignette Quant_13B_problem

sink("Quant_13B_problem_key.txt")
for (i in 1:n) cat(key[i])
cat("\n")
sink()

#. The choice data are input from the .txt file
#. named Quant_13B_problem_U.rda

key <- as.matrix(read.table(file="Quant_13B_problem_key.txt"))

# ----------  Set up 5 choices for each item -------------------------
#. There are only four choices, but a choice for missing or 
#  illegitimate responses is also added.

noption <- matrix(5,n,1)
scrrng  <- c(0,37)

# ----------  print summary count of choices for each option  -------

for (i in 1:n) {
  print(paste("Item: ",i,sep=""))
  for (m in 1:noption[i]) {
    print(paste("Option: ",m," N= ",sum(U[,i]==m),sep=""))
  }
}

# ------ Define the option pre-assigned score values for each item ---

ScoreList <- vector("list",n) # option scores
for (item in 1:n) {
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  ScoreList[[item]] <- scorei
}

optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

# ----------------  Externally define WfdPar object  -----------

Wnbasis <- 2
Wnorder <- min(Wnbasis, 5)
Wbasis  <- fda::create.bspline.basis(c(0,100), Wnbasis, Wnorder) 
WfdPar  <- fdPar(Wbasis)

# ----------------  Initialization Steps  ------------------------

Quant_13B_problem_dataList <- 
  make.dataList(U, key, optList, scrrng, titlestr, NumBasis=2, WfdPar=NULL)

#  save this list object in the data folder

save(Quant_13B_problem_dataList, 
     file="../TestGardener/data/Quant_13B_problem_dataList.rda")  
load(file="../TestGardener/data/Quant_13B_problem_dataList.rda")

#  plot the sum scores as a histogram 

hist(Quant_13B_problem_dataList$scrvec, Quant_13B_problem_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

#  compute the initial option surprisal curves using the 
#  percentage ranks as initial estimates of theta.  This can be left out.

theta     <- Quant_13B_problem_dataList$percntrnk
thetaQnt  <- Quant_13B_problem_dataList$thetaQnt
N         <- length(theta)
n         <- ncol(Quant_13B_problem_dataList$U)

#  preliminary displays and plots before full analysis

H <- Hfun(theta, Quant_13B_problem_dataList$WfdList, U)
meanH <- mean(H)
print(meanH)

WfdResult <- Wbinsmth(theta, Quant_13B_problem_dataList)

#  Plot the initial option proability and surprisal curves

WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)

indfine   <- seq(0,100,len=101)
plotindex <- 1
plotrange <- c(0,100)

ICC.plot(indfine, WfdList, Quant_13B_problem_dataList, Qvec, binctr,  
         Wrng=c(0,3), plotindex=plotindex)

#  ----------------------------------------------------------------------------
#                Proceed through the cycles, returning system time
#  ----------------------------------------------------------------------------

#  Set number of cycles and the cell array to containing the parameter

ncycle <- 10
system.time(
AnalyzeResult <- Analyze(theta, thetaQnt, Quant_13B_problem_dataList, ncycle, 
                         itdisp=TRUE, verbose=TRUE) 
)

#  10 cycles:  mean H = 16.4,  arc length = 53.5
#  system.time = 14.2 seconds

#  set up mean H and arc length for plotting over cycles

HALsave <- matrix(0,ncycle,2)
for (icycle in 1:ncycle) {
  HALsave[icycle,1] <- AnalyzeResult[[icycle]]$meanH
  HALsave[icycle,2] <- AnalyzeResult[[icycle]]$arclength
}

#  ----------------------------------------------------------------------------
#              Plot the average H value and arc length over cycles
#  ----------------------------------------------------------------------------

par(mfrow=c(2,1))

plot(1:ncycle, HALsave[,1], type="b", lwd=2, 
     xlab="Cycle Number",ylab="Mean H")
plot(1:ncycle, HALsave[,2], type="b", lwd=2, 
     xlab="Cycle Number", ylab="Arc Length")

#  select cycle for plotting

icycle <- 10
Quant_13B_problem_parList  <- AnalyzeResult$parList[[icycle]]

#  save this list object in the data folder 5 May 2023

save(Quant_13B_problem_parList, 
     file="../TestGardener/data/Quant_13B_problem_parList.rda") 
load(file="../TestGardener/data/Quant_13B_problem_parList.rda") 

WfdList    <- Quant_13B_problem_parList$WfdList
Qvec       <- Quant_13B_problem_parList$Qvec
binctr     <- Quant_13B_problem_parList$binctr
theta      <- Quant_13B_problem_parList$theta
Hval       <- Quant_13B_problem_parList$Hval
DHval      <- Quant_13B_problem_parList$DHval
D2Hval     <- Quant_13B_problem_parList$D2Hval

#  ----------------------------------------------------------------------------
#                   Compute the arc length or information measure 
#  ----------------------------------------------------------------------------

# A variety of useful objects related to the test info manifold
# are computed and returned in a struct object infoStr

Quant_13B_problem_infoList <- theta2arclen(theta, Qvec, WfdList, binctr)  

save(Quant_13B_problem_infoList, 
     file="../TestGardener/data/Quant_13B_problem_infoList.rda")  
load(file="../TestGardener/data/Quant_13B_problem_infoList.rda")

#  The length of the test manifold
arclength     <- Quant_13B_problem_infoList$arclength
#  indefinite integral of arc length values corresponding to 
#  equally spaced theta values
arclengthvec  <- Quant_13B_problem_infoList$arclengthvec
#  The N arc length values corresponding to the N estimated score indes
#  theta values
theta_al      <- Quant_13B_problem_infoList$theta_al
#  The arc length values for the five marker percentages
Qvec_al       <- Quant_13B_problem_infoList$Qvec_al
#  The arc length values for the bin centers
binctr_al     <- Quant_13B_problem_infoList$binctr_al

#  ----------------------------------------------------------------------------
#        Plot surprisal curves for each test question using ICC.plot
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

indfine   <- seq(0,100,len=101)
plotindex <- 1:n
plotrange <- c(0,100)

#  over score index theta

plotindex=6
ICC.plot(indfine, WfdList, Quant_13B_problem_dataList, Qvec, binctr, 
         data_point=TRUE, plotType=c("P","W"), Wrng=c(0,4), 
         DWrng=c(-0.2, 0.2), plotindex=plotindex)

#  over arclength or information

ICC.plot(arclengthvec, WfdList, Quant_13B_problem_dataList, Qvec_al, binctr_al,  
         data_point=TRUE, plotType=c("P","W"), Wrng=c(0,4))

#  ----------------------------------------------------------------------------
#               Plot density of theta  using scoreDensity
#  ----------------------------------------------------------------------------

ttllab   <- paste(titlestr,": percent rank", sep="")
scrrng   <- c(0,100)
indden10 <- scoreDensity(theta, scrrng, ttlstr=ttllab)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu     <- testscore(theta, WfdList, Quant_13B_problem_dataList$optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- scoreDensity(mu, Quant_13B_problem_dataList$scrrng, ttlstr=ttllab) 

#  compute expected score for each value in the fine mesh of theta values

indfine <- seq(0,100,len=101)
mufine  <- testscore(indfine, WfdList, Quant_13B_problem_dataList$optList)
mu.plot(mufine, Quant_13B_problem_dataList$scrrng, titlestr)

#  ----------------------------------------------------------------------------
#         Plot arc length over a fine mesh of theta values and plot
#  ----------------------------------------------------------------------------

#  print length of the test effort curve

print(paste("Arc length =", round(arclength,2)))

#  plot arc length over fine mesh

ArcLength.plot(arclength, arclengthvec, titlestr)

#  ----------------------------------------------------------------------------
#         plot the distribution of score index values theta 
#  ----------------------------------------------------------------------------

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

par(mfrow=c(2,1))
density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
             titlestr="SweSAT 13B Theta Density",  
             scrnbasis=11, nfine=101)

#  Plot the probability density function as a smooth curve along with
#  points indicating boundary proportions

density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
             titlestr="SweSAT 13B Info Density",  
             scrnbasis=15, nfine=101)

#  ----------------------------------------------------------------------------
#  print test info curve from its first two or three principal components
#  ----------------------------------------------------------------------------

# nharm=2
Results <- Wpca.plot(arclength, WfdList, Quant_13B_problem_dataList$Wdim, 
                     titlestr=titlestr)
print(paste("Variance percent = ",100*round(Results[[3]],3)))

# nharm=3
Results <-Wpca.plot(arclength, WfdList, Quant_13B_problem_dataList$Wdim,   
                    3, dodge = 1.005, titlestr=titlestr)
print(paste("Variance percent = ",100*round(Results[[3]],3)))

#  ----------------------------------------------------------------------------
#                          print sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

Sensitivity.plot(indfine, WfdList, Qvec, Quant_13B_problem_dataList)
  
#  ----------------------------------------------------------------------------
#                          print power curves
#  ----------------------------------------------------------------------------

Power.plot(indfine, WfdList, Qvec, Quant_13B_problem_dataList, height=0.25)

#  ----------------------------------------------------------------------------
#                          print entropy curves
#  ----------------------------------------------------------------------------

Entropy.plot(indfine, WfdList, Qvec, Quant_13B_problem_dataList, height=1)

#  ----------------------------------------------------------------------------
#                   printlay H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

index <- 1:5

Hfuns.plot(indfine, theta, WfdList, Quant_13B_problem_dataList$U, plotindex=100:200)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

simList <- dataSimulation(Quant_13B_problem_dataList, Quant_13B_problem_parList, nsample=50)

scorePerformance(simList, Qvec)

##  Item analysis plots over peak intervals

peakbdry = C(0,21,47,67,100)

Wrng = C(0,5)

#  ----------------------------------------------------------------------------
#                Compute item arc lengths over [0,100]  
#  ----------------------------------------------------------------------------

itemArclength_0_100  <- matrix(0,80,1)
indfine <- seq(0,100,len=1000)
thetai <- NULL
for (i in 1:n) {
  print(i)
  Resulti    <- theta2arclen(theta, Qvec, WfdList, binctr, i)
  thetai     <- Resulti$theta_al
  arclengthi <- Resulti$arclength
  itemArclength_0_100[i] <- arclengthi
}

#  print all of the arc lengths

print("arclength_0_100:")
print(itemArclength_0_100)

#  histogram of the 80 arc lengths 

hist(itemArclength_0_100)

#  sort 80 arclengths in descending order

Result <- sort.int(itemArclength_0_100, decreasing=TRUE, index.return=TRUE)
alsort <- Result$x 
isort  <- Result$ix

#  printlay the four items with the longest arc lengths:
#  28, 63, 24, 26

#  plot the surprisal curves in descending order of arc length

Wrng <- c(0,5)

ICC.plot(WfdList, Quant_13B_problem_dataList, Qvec, binctr, plotType = c("P", "W"), 
         plotindex=isort)

#  plot the surprisal curves in ascending order of arc length

Result <- sort.int(itemArclength_0_100)
alsort <- Result$alsort
isort  <- Resuslt$isort
ICC.plot(WfdList, Quant_13B_problem_dataList, Qvec, binctr, plotType = c("P", "W"), 
         plotindex=isort)

#  ----------------------------------------------------------------------------
#               Compute item arc lengths over the peak intervals 
#  ----------------------------------------------------------------------------

itemArclength_80_100 <- matrix(0,80,1)

indfine <- seq(80,100,len=101)

pltrng = c(80,100)
thetai <- NULL
for (i in 1:80) {
  Result     <- theta2arclen(theta, Qvec, WfdList, binctr, i, pltrng, TRUE)
  thetai     <- Result$theta_al
  arclengthi <- Result$arclength
  itemArclength_80_100[i] <- arclengthi
}

# print('arclength_80_100:')
# print([c((1:80),itemArclength_80_100])

hist(itemArclength_80_100)

Result = sort.int(itemArclength_80_100, decreasing=TRUE, index.return=TRUE)
alsort = Result$x 
isort  = Result$ix

ICC.plot(indfine, WfdList, Quant_13B_problem_dataList, Qvec, binctr, plotType = c("P", "W"), 
         plotindex=isort, plotrange=c(80,100))


