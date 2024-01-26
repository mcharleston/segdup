dat <- read.csv("results.csv", header=T)

d <- 10
l <- 1

par(mfrow=c(1,2))

#compare to Multrec
with(dat, plot(cost ~ mrCost))
abline(0,1,col="red")
with(dat, sum(cost > mrCost))
with(dat, sum(cost < mrCost))
with(dat, sum(dups > mrDups))
with(dat, sum(dups < mrDups))

#compare to true reconciliation
plot(d*dat$nAllDupEvents+l*dat$nLineageSort, dat$cost)
abline(0,1,col="red")
