#Mostly Taken froim figure plotting functions / in association with them

t1 <- table(tnD$v$Time)
t2 <- table(stD$v$Time)
t3 <- table(naD$v$Time)

t4 <- table(tnD_met$v$Time)
t5 <- table(tnD_NM$v$Time)

mean(tnD$e$Distance[tnD$e$Distance < 0.05])
mean(stD$e$Distance[stD$e$Distance < 0.05])
mean(naD$e$Distance[naD$e$Distance < 0.05])

kruskal.test(list(tnD$e$Distance[tnD$e$Distance < 0.05],
                  stD$e$Distance[stD$e$Distance < 0.05],
                  naD$e$Distance[naD$e$Distance < 0.05]))

ks.test(tnD_met$e$Distance[tnD_met$e$Distance < 0.05],
        tnD_NM$e$Distance[tnD_NM$e$Distance < 0.05])

mean(tnD_NM$e$Distance[tnD_NM$e$Distance < 0.05])
mean(tnD_met$e$Distance[tnD_met$e$Distance < 0.05])

kruskal.test(list(tnD$e$Distance,
                  stD$e$Distance,
                  naD$e$Distance))

ks.test(tnD_met$e$Distance,
        tnD_NM$e$Distance)                  

mean(tnD$e$Distance)
mean(stD$e$Distance)
mean(naD$e$Distance)

mean(tnD_met$e$Distance)
mean(tnD_NM$e$Distance)


ageD <- lapply(tnD, function(x) x$f)
ageDi <- ageD[["0.04"]]
mod <- glm(cbind(Positive, vTotal) ~ tDiff+ oeDens, data=ageDi, family='binomial')
1-summary(mod)$deviance/summary(mod)$null.deviance

ageD <- lapply(stD, function(x) x$f)
ageDi <- ageD[["0.04"]]
mod <- glm(cbind(Positive, vTotal) ~ tDiff+ oeDens, data=ageDi, family='binomial')
1-summary(mod)$deviance/summary(mod)$null.deviance

ageD <- lapply(naD, function(x) x$f)
ageDi <- ageD[["0.04"]]
mod <- glm(cbind(Positive, vTotal) ~ tDiff+ oeDens, data=ageDi, family='binomial')
1-summary(mod)$deviance/summary(mod)$null.deviance

ageD <- lapply(tnD_NM, function(x) x$f)
ageDi <- ageD[["0.04"]]
mod <- glm(cbind(Positive, vTotal) ~ tDiff+ oeDens, data=ageDi, family='binomial')
1-summary(mod)$deviance/summary(mod)$null.deviance

ageD <- lapply(tnD_met, function(x) x$f)
ageDi <- ageD[["0.04"]]
ageDi <- subset(ageDi, vTotal>63 & tDiff<14)
mod <- glm(cbind(Positive, vTotal) ~ tDiff + oeDens, data=ageDi, family='binomial')
1-summary(mod)$deviance/summary(mod)$null.deviance

dMax <- max(subG$e$Distance)
vTab <- table(subG$v$Time)
eTab <- rep(0, length(vTab))
names(eTab) <- names(vTab)
temp1 <- table(subG$e$t1) 
temp2 <- table(subG$e$t2)
temp3 <- table(subset(subG$e, t1==t2)$tMax)
eTab[names(temp1)] <- eTab[names(temp1)]+unname(temp1)
eTab[names(temp2)] <- eTab[names(temp2)]+unname(temp2)
eTab[names(temp3)] <- eTab[names(temp3)]-unname(temp3)
plot(eTab/vTab)


ageD <- lapply(tnD, function(x) x$g)


runs <- readRDS("~/Data/Paper1/tn93StsubB_RD.rds")
r3 <- readRDS("~/Data/Paper1/tn93TnsubB_nomet_RD.rds")
r4 <- readRDS("~/Data/Paper1/tn93TnsubB_met_RD.rds")

