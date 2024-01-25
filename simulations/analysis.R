dat <- read.csv("results.csv", header=T)

d <- 10
l <- 1

with(dat, plot(cost ~ mrCost, ylim=c(0,200)))
abline(0,1,col="red")

plot(d*dat$nAllDupEvents+l*dat$nLineageSort, dat$cost, ylim=c(0,200))
abline(0,1,col="red")
