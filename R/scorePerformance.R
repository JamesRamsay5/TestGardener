scorePerformance <- function(dataList, simList) {

sumscrsave <- simList$sumscrsave
indexsave  <- simList$indexsave
musave     <- simList$musave
infosave   <- simList$infosave
index.pop  <- simList$index.pop
mu.pop     <- simList$mu.pop
info.pop   <- simList$info.pop

nindex.pop <- dim(indexsave)[1]
nsample    <- dim(indexsave)[2]

indexres  <- indexsave  - matrix(index.pop, nindex.pop, nsample)
mures     <- musave     - matrix(mu.pop,    nindex.pop, nsample)
sumscrres <- sumscrsave - matrix(mu.pop,    nindex.pop, nsample)

#  compute bias's averaged over only samples

indexbias  <- apply(indexres,  1, mean)
mubias     <- apply(mures,     1, mean)
sumscrbias <- apply(sumscrres, 1, mean)

#  compute RMSE's averaged over only samples

indexRMSE  <- sqrt(apply( indexres^2, 1, mean))
muRMSE     <- sqrt(apply(    mures^2, 1, mean))
sumscrRMSE <- sqrt(apply(sumscrres^2, 1, mean))

#  compute Variance's averaged over only samples

indexstd  <- sqrt( indexRMSE^2  -  indexbias^2)
mustd     <- sqrt(    muRMSE^2  -     mubias^2)
sumscrstd <- sqrt(sumscrRMSE^2  - sumscrbias^2)

#  Set up bases for returning smooth performance curves

indbasis <- create.bspline.basis(c(0,100), 13)  
indfdPar <- fdPar(indbasis, 2, 1e0)

#  Lightly smooth bias results over index's

indexbiasfd  <- smooth.basis(index.pop,     indexbias, indfdPar)$fd
mubiasfd     <- smooth.basis(index.pop,        mubias, indfdPar)$fd
sumscrbiasfd <- smooth.basis(index.pop,    sumscrbias, indfdPar)$fd

#  Lightly smooth RMSE results over index's

indexRMSEfd  <- smooth.basis(index.pop,     indexRMSE, indfdPar)$fd
muRMSEfd     <- smooth.basis(index.pop,        muRMSE, indfdPar)$fd
sumscrRMSEfd <- smooth.basis(index.pop,    sumscrRMSE, indfdPar)$fd
#  Lightly smooth std. dev. results over index's

indexstdfd   <- smooth.basis(index.pop,     indexstd, indfdPar)$fd
mustdfd      <- smooth.basis(index.pop,        mustd, indfdPar)$fd
sumscrstdfd  <- smooth.basis(index.pop,    sumscrstd, indfdPar)$fd

#  get results for arc length if required

if (!is.null(infosave)) {
    infores    <- infosave - matrix(info.pop,nindex.pop,nsample)
    infobias   <- apply(infores, 1, mean)
    infoRMSE   <- sqrt(apply(infores^2, 1, mean))
    infostd    <- sqrt(infoRMSE^2 - infobias^2)
    infobiasfd <- smooth.basis(index.pop, infobias, indfdPar)$fd
    infoRMSEfd <- smooth.basis(index.pop, infoRMSE, indfdPar)$fd
    infostdfd  <- smooth.basis(index.pop, infostd,  indfdPar)$fd
} else {
    infobiasfd <- NULL
    infoRMSEfd <- NULL
    infostdfd  <- NULL
}

#  bundle these results into struct object simList

simListout <- simList

indfine <- seq(0,100,len=101)

simListout$indexRMSE  <- eval.fd(indfine, indexRMSEfd)
simListout$muRMSE     <- eval.fd(indfine, muRMSEfd)
simListout$sumscrRMSE <- eval.fd(indfine, sumscrRMSEfd)
simListout$indexbias  <- eval.fd(indfine, indexbiasfd)
simListout$mubias     <- eval.fd(indfine, mubiasfd)
simListout$sumscrbias <- eval.fd(indfine, sumscrbiasfd)
simListout$indexstd   <- eval.fd(indfine, indexstdfd)
simListout$mustd      <- eval.fd(indfine, mustdfd)
simListout$sumscrstd  <- eval.fd(indfine, sumscrstdfd)

if (!is.null(infosave)) {
    simListout$infoRMSE <- eval.fd(indfine, infoRMSEfd)
    simListout$infobias <- eval.fd(indfine, infobiasfd)
    simListout$infostd  <- eval.fd(indfine, infostdfd)
} else {
    simListout$infoRMSE <- NULL
    simListout$infobias <- NULL
    simListout$infostd  <- NULL
}

 return(simListout)

}
