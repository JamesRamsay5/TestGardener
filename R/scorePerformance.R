scorePerformance <- function(dataList, simList) {

sumscrsave <- simList$sumscr
thetasave  <- simList$theta
musave     <- simList$mu
alsave     <- simList$al
theta.pop  <- simList$thepop
mu.pop     <- simList$mupop
al.pop     <- simList$alpop

ntheta.pop <- dim(thetasave)[1]
nsample    <- dim(thetasave)[2]

thetares  <- thetasave  - matrix(theta.pop,ntheta.pop,nsample)
mures     <- musave     - matrix(mu.pop,   ntheta.pop,nsample)
sumscrres <- sumscrsave - matrix(mu.pop,   ntheta.pop,nsample)

#  compute bias's averaged over only samples

thetabias  <- apply(thetares,  1, mean)
mubias     <- apply(mures,     1, mean)
sumscrbias <- apply(sumscrres, 1, mean)

#  compute RMSE's averaged over only samples

thetaRMSE  <- sqrt(apply( thetares^2, 1, mean))
muRMSE     <- sqrt(apply(    mures^2, 1, mean))
sumscrRMSE <- sqrt(apply(sumscrres^2, 1, mean))

#  compute Variance's averaged over only samples

thetastd  <- sqrt( thetaRMSE^2  -  thetabias^2)
mustd     <- sqrt(    muRMSE^2  -     mubias^2)
sumscrstd <- sqrt(sumscrRMSE^2  - sumscrbias^2)

#  Set up bases for returning smooth performance curves

indbasis <- create.bspline.basis(c(0,100), 13)  
indfdPar <- fdPar(indbasis, 2, 1e0)

#  Lightly smooth bias results over theta's

thetabiasfd  <- smooth.basis(theta.pop,     thetabias, indfdPar)$fd
mubiasfd     <- smooth.basis(theta.pop,        mubias, indfdPar)$fd
sumscrbiasfd <- smooth.basis(theta.pop,    sumscrbias, indfdPar)$fd

#  Lightly smooth RMSE results over theta's

thetaRMSEfd  <- smooth.basis(theta.pop,     thetaRMSE, indfdPar)$fd
muRMSEfd     <- smooth.basis(theta.pop,        muRMSE, indfdPar)$fd
sumscrRMSEfd <- smooth.basis(theta.pop,    sumscrRMSE, indfdPar)$fd
#  Lightly smooth std. dev. results over theta's

thetastdfd   <- smooth.basis(theta.pop,     thetastd, indfdPar)$fd
mustdfd      <- smooth.basis(theta.pop,        mustd, indfdPar)$fd
sumscrstdfd  <- smooth.basis(theta.pop,    sumscrstd, indfdPar)$fd

#  get results for arc length if required

if (!is.null(alsave)) {
    alres    <- alsave - matrix(al.pop,ntheta.pop,nsample)
    albias   <- apply(alres, 1, mean)
    alRMSE   <- sqrt(apply(alres^2, 1, mean))
    alstd    <- sqrt(alRMSE^2 - albias^2)
    albiasfd <- smooth.basis(theta.pop, albias, indfdPar)$fd
    alRMSEfd <- smooth.basis(theta.pop, alRMSE, indfdPar)$fd
    alstdfd  <- smooth.basis(theta.pop, alstd,  indfdPar)$fd
} else {
    albiasfd <- NULL
    alRMSEfd <- NULL
    alstdfd  <- NULL
}

#  bundle these results into struct object simList

simListout <- simList

indfine <- seq(0,100,len=101)

simListout$thetaRMSE  <- eval.fd(indfine, thetaRMSEfd)
simListout$muRMSE     <- eval.fd(indfine, muRMSEfd)
simListout$sumscrRMSE <- eval.fd(indfine, sumscrRMSEfd)
simListout$thetabias  <- eval.fd(indfine, thetabiasfd)
simListout$mubias     <- eval.fd(indfine, mubiasfd)
simListout$sumscrbias <- eval.fd(indfine, sumscrbiasfd)
simListout$thetastd   <- eval.fd(indfine, thetastdfd)
simListout$mustd      <- eval.fd(indfine, mustdfd)
simListout$sumscrstd  <- eval.fd(indfine, sumscrstdfd)

if (!is.null(alsave)) {
    simListout$alRMSE <- eval.fd(indfine, alRMSEfd)
    simListout$albias <- eval.fd(indfine, albiasfd)
    simListout$alstd  <- eval.fd(indfine, alstdfd)
} else {
    simListout$alRMSE <- NULL
    simListout$albias <- NULL
    simListout$alstd  <- NULL
}

 return(simListout)

}
