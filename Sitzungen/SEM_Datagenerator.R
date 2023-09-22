library(lavaan)
source("SEMgraph.R")

# simulation model definition ----
simModel = '
# 1. Faktoren

sleep1 =~ psqi1 + 0.5*sh1 + 0.8*fas1
sleep2 =~ psqi2 + 0.5*sh2 + 0.8*fas2
sleep3 =~ psqi3 + 0.5*sh3 + 0.8*fas3
sleep4 =~ psqi4 + 0.5*sh4 + 0.8*fas4
mwb1 =~ q1mwb1 + 0.7*q2mwb1 + 0.4*q3mwb1
mwb2 =~ q1mwb2 + 0.7*q2mwb2 + 0.4*q3mwb2
mwb3 =~ q1mwb3 + 0.7*q2mwb3 + 0.4*q3mwb3
mwb4 =~ q1mwb4 + 0.7*q2mwb4 + 0.4*q3mwb4
isleep =~ 1 * sleep1 + 1 * sleep2 + 1 * sleep3 + 1 * sleep4
imwb =~ 1 * mwb1 + 1 * mwb2 + 1 * mwb3 + 1 * mwb4

# 2. Autoregression
sleep2 ~ 0.3 * sleep1
sleep3 ~ 0.3 * sleep2
sleep4 ~ 0.3 * sleep3
mwb2 ~ 0.1 * mwb1
mwb3 ~ 0.1 * mwb2
mwb4 ~ 0.1 * mwb3

# 3. Kreuzregression
mwb2 ~ 0.11 * sleep1
mwb3 ~ 0.11 * sleep2
mwb4 ~ 0.11 * sleep3
sleep2 ~ 0.46 * mwb1
sleep3 ~ 0.46 * mwb2
sleep4 ~ 0.46 * mwb3

# 4. Kovarianzen
sleep1 ~~ 0.1 * mwb1
sleep2 ~~ 0.1 * mwb2
sleep3 ~~ 0.1 * mwb3
sleep4 ~~ 0.1 * mwb4

isleep ~~ 0.6 * imwb

sleep1 ~~ 0 * isleep
mwb1 ~~ 0 * isleep
sleep1 ~~ 0 * imwb
mwb1 ~~ 0 * imwb

# 5. Intercepts
isleep ~ 1.08 * 1
imwb ~ 0.67 * 1
'
# ----

dat = simulateData(simModel, sample.nobs=2837)

rescale = function(vec, newMean, newSD) { scale(vec) * newSD + newMean }
cols = c("sh1", "sh2", "sh3", "sh4")
dat[cols] = lapply(dat[cols], rescale, newMean=6.85, newSD=0.64)

# clpmri definition ----
clpmri = '
sleep1 =~ psqi1 + sl2 * sh1 + sl3 * fas1
sleep2 =~ psqi2 + sl2 * sh2 + sl3 * fas2
sleep3 =~ psqi3 + sl2 * sh3 + sl3 * fas3
sleep4 =~ psqi4 + sl2 * sh4 + sl3 * fas4
mwb1 =~ q1mwb1 + mwbl2 * q2mwb1 + mwbl3 * q3mwb1
mwb2 =~ q1mwb2 + mwbl2 * q2mwb2 + mwbl3 * q3mwb2
mwb3 =~ q1mwb3 + mwbl2 * q2mwb3 + mwbl3 * q3mwb3
mwb4 =~ q1mwb4 + mwbl2 * q2mwb4 + mwbl3 * q3mwb4
isleep =~ 1 * sleep1 + 1 * sleep2 + 1 * sleep3 + 1 * sleep4
imwb =~ 1 * mwb1 + 1 * mwb2 + 1 * mwb3 + 1 * mwb4


# 2. Autoregression
sleep2 ~ a1 * sleep1
sleep3 ~ a1 * sleep2
sleep4 ~ a1 * sleep3
mwb2 ~ a2 * mwb1
mwb3 ~ a2 * mwb2
mwb4 ~ a2 * mwb3

# 3. Kreuzregression
mwb2 ~ c1 * sleep1
mwb3 ~ c1 * sleep2
mwb4 ~ c1 * sleep3
sleep2 ~ c2 * mwb1
sleep3 ~ c2 * mwb2
sleep4 ~ c2 * mwb3

# 4. Kovarianzen
sleep1 ~~ cov1 * mwb1
sleep2 ~~ cov2 * mwb2
sleep3 ~~ cov2 * mwb3
sleep4 ~~ cov2 * mwb4

isleep ~~ covi*imwb

sleep1 ~~ 0 * isleep
mwb1 ~~ 0 * isleep
sleep1 ~~ 0 * imwb
mwb1 ~~ 0 * imwb

# 5. Intercepts
psqi1 ~ 0 * 1
psqi2 ~ 0 * 1
psqi3 ~ 0 * 1
psqi4 ~ 0 * 1
sh1 ~ 1
sh2 ~ 1
sh3 ~ 1
sh4 ~ 1
fas1 ~ 0 * 1
fas2 ~ 0 * 1
fas3 ~ 0 * 1
fas4 ~ 0 * 1
q1mwb1 ~ 0 * 1
q1mwb2 ~ 0 * 1
q1mwb3 ~ 0 * 1
q1mwb4 ~ 0 * 1
q2mwb1 ~ 0 * 1
q2mwb2 ~ 0 * 1
q2mwb3 ~ 0 * 1
q2mwb4 ~ 0 * 1
q3mwb1 ~ 0 * 1
q3mwb2 ~ 0 * 1
q3mwb3 ~ 0 * 1
q3mwb4 ~ 0 * 1

isleep ~ 1
imwb ~ 1
'
#----
clpmrifit = sem(clpmri, data=dat)
summary(clpmrifit, fit.measures=T)
SEMgraph(clpmrifit) # for funzies

# If the data looks good, save it as csv to load it again in the next section
# write.csv(dat, "SleepMWB.csv", row.names=F)

#----------------------------------------------------------------------------

#library(psych)
rm(dat, clpmrifit)
dat = read.csv("SleepMWB.csv")
clpmrifit = sem(clpmri, data=dat)
factorLoadings = parameterEstimates(clpmrifit)
factorLoadings = factorLoadings[factorLoadings$op == "=~" & !(factorLoadings$lhs %in% c("isleep", "imwb")),]$est
loadingMatrix = matrix(0, nrow=24, ncol=8)
for(i in 1:(length(factorLoadings) / 3)) {
  startRow = (i-1)*3+1
  loadingMatrix[startRow:(startRow+2), i] = factorLoadings[startRow:(startRow+2)]
}

#scoreDat = data.frame(factor.scores(dat, loadingMatrix, method="tenBerge")$scores) # needs psych package
#scoreDat = data.frame(as.matrix(dat) %*% (cor(dat) %*% loadingMatrix))
scoreDat = as.data.frame(lavPredict(clpmrifit)[,1:8]) # most appropriate technique
colnames(scoreDat) = c("fssleep1", "fssleep2", "fssleep3", "fssleep4", "fsmwb1", "fsmwb2", "fsmwb3", "fsmwb4")

# clpmrishort definition ----
clpmrishort = '
sleep1 =~ 1 * fssleep1
sleep2 =~ 1 * fssleep2
sleep3 =~ 1 * fssleep3
sleep4 =~ 1 * fssleep4
mwb1 =~ 1 * fsmwb1
mwb2 =~ 1 * fsmwb2
mwb3 =~ 1 * fsmwb3
mwb4 =~ 1 * fsmwb4
isleep =~ 1 * sleep1 + 1 * sleep2 + 1 * sleep3 + 1 * sleep4
imwb =~ 1 * mwb1 + 1 * mwb2 + 1 * mwb3 + 1 * mwb4


# 2. Autoregression
sleep2 ~ a1 * sleep1
sleep3 ~ a1 * sleep2
sleep4 ~ a1 * sleep3
mwb2 ~ a2 * mwb1
mwb3 ~ a2 * mwb2
mwb4 ~ a2 * mwb3

# 3. Kreuzregression
mwb2 ~ c1 * sleep1
mwb3 ~ c1 * sleep2
mwb4 ~ c1 * sleep3
sleep2 ~ c2 * mwb1
sleep3 ~ c2 * mwb2
sleep4 ~ c2 * mwb3

# 4. Kovarianzen
sleep1 ~~ cov1 * mwb1
sleep2 ~~ cov2 * mwb2
sleep3 ~~ cov2 * mwb3
sleep4 ~~ cov2 * mwb4

isleep ~~ covi*imwb

sleep1 ~~ 0 * isleep
mwb1 ~~ 0 * isleep
sleep1 ~~ 0 * imwb
mwb1 ~~ 0 * imwb

# 5. Intercepts
fssleep1 ~ 0 * 1
fssleep2 ~ 0 * 1
fssleep3 ~ 0 * 1
fssleep4 ~ 0 * 1
fsmwb1 ~ 0 * 1
fsmwb2 ~ 0 * 1
fsmwb3 ~ 0 * 1
fsmwb4 ~ 0 * 1

isleep ~ 1
imwb ~ 1
'
#----

clpmrishortfit = sem(clpmrishort, scoreDat)
summary(clpmrishortfit, fit.measures=T)
SEMgraph(clpmrishortfit)
pebig = parameterEstimates(clpmrifit)
pesmall = parameterEstimates(clpmrishortfit)
merged = merge(pebig[c(cols, "est")], pesmall[c(cols, "est")], by=cols)
estimateError = sqrt(sum((merged$est.x - merged$est.y)^2))

#write.csv(scoreDat, "SleepMWBScored.csv", row.names=F)