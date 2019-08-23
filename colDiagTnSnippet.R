runsMet <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_met_LD.rds")
runsNM <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_nomet_LD.rds")
table(runsMet[[5]]$`0`$v$Time)
table(runsNM[[5]]$`0`$v$Time)

runsMet <- rev(runsMet)
runsNM <- rev(runsNM)

#Obtain a list of vectors of GAICs for each filtered run
cutoffs <- as.numeric(names(runsNM[[1]]))
#The step distance between cutoff points
step <- max(cutoffs) / (length(cutoffs)-1)

gaicsMet <- lapply(runsMet, function(run){sapply(run, function(x) {x$gaic})})

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLocMet <- sapply(gaicsMet, function(x){step*(which(x==min(x))[[1]]-1)}) 
minsMet <- sapply(gaicsMet, function(x){min(x)}) 
minsMet

gaicsNM <- lapply(runsNM, function(run){sapply(run, function(x) {x$gaic})})

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLocNM <- sapply(gaicsNM, function(x){step*(which(x==min(x))[[1]]-1)}) 
minsNM <- sapply(gaicsNM, function(x){min(x)}) 
minsNM

yMet <- tail(unname(table(runsMet[[5]]$`0`$v$Time)), 5)
yNM <- tail(unname(table(runsNM[[5]]$`0`$v$Time)), 5)

plot(c(yMet, yNM),c(minsMet,minsNM)) 
