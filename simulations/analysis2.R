dat <- read.csv("allresults.csv", header=T)

#default values (change nP to 20)
nH <- 50
nP <- 30
rB <- 2
pJ <- 0.5

d <- 10
l <- 1
dat$trueCost <- d*dat$nAllDupEvents+l*dat$nLineageSort

dat$propSdCost <- (dat$cost - dat$mrCost)/dat$mrCost
dat$propTrueCost <- (dat$trueCost - dat$cost)/dat$cost

dat$sdBetter <- dat$cost < dat$mrCost
dat$same <- dat$cost == dat$mrCost
dat$mrBetter <- dat$mrCost < dat$cost

#varying nH
nHdat <- dat[dat$nP==nP & dat$rB==rB & dat$pJ == pJ,]

totals <- aggregate(nHdat, by = list(nHdat$nH), FUN = sum)
means <- aggregate(nHdat, by = list(nHdat$nH), FUN = mean)
sds <- aggregate(nHdat, by = list(nHdat$nH), FUN = sd)

#compare cost to Multrec - proportional increase
plot(means$nH, means$propSdCost)

#compare cost to Multrec - counts
barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$nH, beside=T, col=c("red","green","blue"))

#compare time to Multrec
plot(means$nH, means$sdTime, col="red")
points(means$nH, means$mrTime, col="blue")

#compare cost to true reconciliation
plot(means$nH, means$propTrueCost)
