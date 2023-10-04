library(TestGardener)

# build and install TestGardener from the revised code here

setwd("../UManitoba/hads")

library(plotly)
library(dplyr)
library(ggplot2)

#  ----------------------------------------------------------------------------
##   load versions of dataList, parList and infoList
#  ----------------------------------------------------------------------------

#. load n, U and dataList for setting up for analysis

n <- 7
load(file="U_dep_3.rda")
load(file="hads_dep_3d_4b_dataList.rda")

#. load parList and infoList 

load(file="hads_dep_3d_4b_parList.rda")
load(file="hads_dep_3d_4b_infoList.rda")

#  ----------------------------------------------------------------------------
#         Select the set of parList objects
#  ----------------------------------------------------------------------------

#. four basis functions

WfdList_4b    <- hads_dep_3d_4b_parList$WfdList
meanH_4b      <- hads_dep_3d_4b_parList$meanH
arclength_4b  <- hads_dep_3d_4b_parList$arclength
Qvec_4b       <- hads_dep_3d_4b_parList$Qvec
binctr_4b     <- hads_dep_3d_4b_parList$binctr
theta_4b      <- hads_dep_3d_4b_parList$theta

#  ----------------------------------------------------------------------------
#           Select the set of infoList objects 
#  ----------------------------------------------------------------------------

#. four basis functions

arclength_4b     <- hads_dep_3d_4b_infoList$arclength
arclengthvec_4b  <- hads_dep_3d_4b_infoList$arclengthvec
theta_al_4b      <- hads_dep_3d_4b_infoList$theta_al
Qvec_al_4b       <- hads_dep_3d_4b_infoList$Qvec_al
binctr_al_4b     <- hads_dep_3d_4b_infoList$binctr_al

#  ----------------------------------------------------------------------------
#           Plot Info Curve using Wpca.plot.R 
#  ----------------------------------------------------------------------------

titlestr_dep <- "Depression Info"
library(rgl)
options(rgl.useNULL = TRUE) # Suppress the separate window.

Result <- TestGardener::Wpca.plot(arclength_4b, WfdList_4b, 
                                  hads_dep_3d_4b_dataList$Wdim, 
                                  3, titlestr=titlestr_dep, rotate=FALSE)
p <- Result$pcaplt
class(p)
print(p)

#  ----------------------------------------------------------------------------
#                   plot H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

#. over score index

# indfine <- seq(0,100,len=101)
# 
# nindex <- c(1:7)
# 
# evalarg <- indfine
# Nindex <- 1
# theta   <- theta_4b[Nindex]
# WfdList <- WfdList_4b
# U       <- hads_dep_3d_4b_dataList$U[Nindex,]
# 
# Hfine <- Hcurve(WfdList, U)
# plot(evalarg, Hfine, type="l")
# ylim <- c(min(Hfine), max(Hfine))
# lines(c(theta,theta), ylim, type="l", lty=2)
# 
# DHfine <- Hcurve(WfdList, U, nderiv=1)
# plot(evalarg, DHfine, type="l")
# ylim <- c(min(DHfine), max(DHfine))
# lines(c(theta,theta), ylim, type="l", lty=2)
# ylim <- c(0,0)
# lines(c(0,100), ylim, type="l", lty=2)
# 
# D2Hfine <- Hcurve(WfdList, U, nderiv=2)
# plot(evalarg, D2Hfine, type="l")
# ylim <- c(min(D2Hfine), max(D2Hfine))
# lines(c(theta,theta), ylim, type="l", lty=2)
# ylim <- c(min(D2Hfine),min(D2Hfine))
# lines(c(0,100), ylim, type="l", lty=2)
# 
# #  ----------------------------------------------------------------------------
# #         plot H and DH curves for raters 1 and 16 in four panels
# #  ----------------------------------------------------------------------------
# 
# par(mfcol=c(2,2))
# 
# #. first rater, H curve
# U1     <- hads_dep_3d_4b_dataList$U[1,]
# theta1 <- theta_4b[1]
# Hfine1 <- Hcurve(WfdList, U1)
# plot(evalarg, Hfine1, type="l", 
#      xlab="", ylab="Fitting Criterion H", ylim=c(0,25))
# ylim <- c(0,25)
# lines(c(theta1,theta1), ylim, type="l", lty=2)
# #. first rater, DH curve
# DHfine1 <- Hcurve(WfdList, U1, nderiv=1)
# plot(evalarg, DHfine1, type="l", 
#      xlab="Score Index", ylab="Fitting Criterion Slope DH", ylim=c(-0.5,1.5))
# ylim <- c(-0.5,1.5)
# lines(c(theta1,theta1), ylim, type="l", lty=2)
# ylim <- c(0,0)
# lines(c(0,100), ylim, type="l", lty=2)
# 
# #. 16th rater, H curve
# U16     <- hads_dep_3d_4b_dataList$U[16,]
# theta16 <- theta_4b[16]
# Hfine16 <- Hcurve(WfdList, U16)
# plot(evalarg, Hfine16, type="l", 
#      xlab="", ylab="Fitting Criterion H", ylim=c(0,25))
# ylim <- c(0,25)
# lines(c(theta16,theta16), ylim, type="l", lty=2)
# #. first rater, DH curve
# DHfine16 <- Hcurve(WfdList, U16, nderiv=1)
# plot(evalarg, DHfine16, type="l", 
#      xlab="Score Index", ylab="Fitting Criterion Slope DH", ylim=c(-0.5,1.5))
# ylim <- c(-0.5,1.5)
# lines(c(theta16,theta16), ylim, type="l", lty=2)
# ylim <- c(0,0)
# lines(c(0,100), ylim, type="l", lty=2)
# 
# #  ----------------------------------------------------------------------------
# #. plot negative log likelihood criterion for first 100 raters
# #  ----------------------------------------------------------------------------
# 
# Nindex <- 1:100
# for (index in Nindex) {
#   U       <- hads_dep_3d_4b_dataList$U[index,]
#   Hfine <- Hcurve(WfdList, U)
#   theta   <- theta_4b[index]
#   plot(evalarg, Hfine, type="l", lwd=4, 
#        xlab="Score Index", ylab="Surprisal (3-bits)",
#        main = paste("Rater",index))
#   ylim <- c(0, max(Hfine))
#   lines(c(theta,theta), ylim, type="l", lty=2, lwd=2)
#   if (length(Nindex) > 1) {
#     readline(prompt = paste("theta", index, ". Press [enter] to continue"))
#   }
# }
# 
# #. plot first derivative of negative log likelihood criterion
# 
# Nindex <- 1:100
# for (index in Nindex) {
#   U       <- hads_dep_3d_4b_dataList$U[index,]
#   DHfine   <- Hcurve(WfdList, U, nderiv=1)
#   theta   <- theta_4b[index]
#   plot(evalarg, DHfine, type="l", lwd=4, 
#        xlab="Score Index", ylab="Surprisal First Derivative (3-bits)",
#        main = paste("Rater",index))
#   ylim <- c(min(DHfine), max(DHfine))
#   lines(c(theta,theta), ylim, type="l", lty=2, lwd=2)
#   ylim <- c(0,0)
#   lines(c(0,100), ylim, type="l", lty=2, lwd=2)
#   if (length(Nindex) > 1) {
#     readline(prompt = paste("theta", index, ". Press [enter] to continue"))
#   }
# }
# 
# #. plot second derivative of negative log likelihood criterion
# 
# Nindex <- 1:100
# for (index in Nindex) {
#   U       <- hads_dep_3d_4b_dataList$U[index,]
#   D2Hfine   <- Hcurve(WfdList, U, nderiv=2)
#   theta   <- theta_4b[index]
#   plot(evalarg, D2Hfine, type="l", lwd=4, 
#        xlab="Score Index", ylab="Surprisal Second Derivative (3-bits)",
#        main = paste("Rater",index))
#   ylim <- c(min(D2Hfine), max(D2Hfine))
#   lines(c(theta,theta), ylim, type="l", lty=2, lwd=2)
#   ylim <- c(min(D2Hfine),min(D2Hfine))
#   lines(c(0,100), ylim=, type="l", lty=2, lwd=2)
#   if (length(Nindex) > 1) {
#     readline(prompt = paste("theta", index, ". Press [enter] to continue"))
#   }
# }
# 
# #. Here is where the failure at the invoking of ggarrange() occurs.
# 
# # error message is:
# 
# # p <- ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)
# # Error in `purrr::map()`:
# #   ℹ In index: 1.
# # Caused by error in `ggplot2::geom_line()`:
# #   ! Problem while computing aesthetics.
# # ℹ Error occurred in the 1st layer.
# # Caused by error in `check_aesthetics()`:
# #   ! Aesthetics must be either length 1 or the same as the data (101)
# # ✖ Fix the following mappings: `y`
# # Run `rlang::last_trace()` to see where the error occurred.
# 
# plotindex = 1
# Hfuns.plot(evalarg, theta, WfdList, U, plotindex)
# 
# evalarg <- arclengthvec_4b
# Nindex  <- 1:5
# theta   <- theta_al_4b[Nindex]
# U       <- hads_dep_3d_4b_dataList$U[Nindex,]
# Hfine   <- Hcurve(evalarg, theta, WfdList, U, nderiv=0)
# 
# matplot(evalarg, Hfine, type="l")
# 
# Hfuns.plot(evalarg,      theta,   WfdList, U, plotindex)
# 
# Hfuns.plot(arclengthvec_4b, theta_al_4b, WfdList_4b, U, 
#            plotindex=1:5)
# plot(evalarg, Hfine, type="l")
# 
# 
