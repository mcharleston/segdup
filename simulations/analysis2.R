dat <- read.csv("allresults.csv", header=T)

nH <- 50
nP <- 20
rB <- 2
pJ <- 0.5

d <- 10
l <- 1

#varying nH
nHdat <- dat[dat$nP==nP & dat$rB==rB & dat$pJ == pJ,]

#compare cost to Multrec
with(dat, plot(cost ~ mrCost))
abline(0,1,col="red")

#compare time to Multrec
with(dat, plot(sdTime ~ mrTime))
abline(0,1,col="red")

#compare cost to true reconciliation
plot(d*dat$nAllDupEvents+l*dat$nLineageSort, dat$cost)
abline(0,1,col="red")
