##TO-DO: Specify source-file Location
source("~/git/tn/comp_An.R")
require(R.utils)

gMet <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_met_G.rds")
gNM <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_G.rds")

tMet <- table(gMet$v$Time)
tNM <- table(gNM$v$Time)

mMet <- as.numeric(max(names(tMet[tMet>100])))
mNM <- as.numeric(max(names(tNM[tNM>100])))

gMet <- tFilt(gMet, mMet)
gNM <- tFilt(gNM, mNM)

tMet <- table(gMet$v$Time)
tNM <- table(gNM$v$Time)


remV <- rev(gNM$v$ID)[1:(153-125)]
gNM$v <- subset(gNM$v, !ID%in%remV)
gNM$e <- subset(gNM$e, !ID1%in%remV & !ID2%in%remV)

gMet <- clsFilt(gMet)
gNM <- clsFilt(gNM)

gMet$f <- likData(gMet)
gNM$f <- likData(gNM)

saveRDS(gMet, "compareMet.rds")
saveRDS(gNM, "compareNM.rds")
