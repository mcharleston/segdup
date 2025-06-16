library(ggplot2)
library(reshape2)
theme_set(theme_bw(base_size=20))

#system("python3 process.py")
dat <- read.csv("allresults-long.csv", header=T)

#default values (change nP to 20)
nH <- 50
nP <- 20
rB <- 2
pJ <- 0.5

d <- 10
l <- 1
dat$trueCost <- d*dat$nAllDupEvents+l*dat$nLineageSort

dat$propSdCost <- dat$cost/dat$mrCost
dat$propTrueCost <- dat$cost/dat$trueCost

dat$sdBetter <- dat$cost < dat$mrCost
dat$same <- dat$cost == dat$mrCost
dat$mrBetter <- dat$mrCost < dat$cost


#varying nH
nHdat <- dat[dat$nP==nP & dat$rB==rB & dat$pJ == pJ,]

totals <- aggregate(nHdat, by = list(nHdat$nH), FUN = sum)
means <- aggregate(nHdat, by = list(nHdat$nH), FUN = mean)
sds <- aggregate(nHdat, by = list(nHdat$nH), FUN = sd)/10

#compare cost to Multrec - proportional increase
#pdf("figures/nH-sdvmr-cost-long.pdf")
#plot(means$nH, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="nH", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
#arrows(means$nH, means$propSdCost-2*sds$propSdCost, means$nH, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
#dev.off()

pdf("figures/nH-sdvmr-cost-long.pdf")
ggplot(data=means) + geom_point(aes(nH,propSdCost)) +
	geom_errorbar(aes(nH,ymin=propSdCost-2*sds$propSdCost,ymax=propSdCost+2*sds$propSdCost)) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab("nH") + ylab("Proportional cost")
dev.off()

#compare cost to Multrec - counts
#pdf("figures/nH-sdvmr-counts-long.pdf")
#barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$nH, beside=T, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="nH")
#legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
#dev.off()

tm <- melt(totals[,c("nH","sdBetter","same","mrBetter")], id.vars=1)
pdf("figures/nH-sdvmr-counts-long.pdf")
ggplot(data=tm) + geom_bar(aes(nH/100,value, fill=variable), position="dodge", stat="identity") +
	scale_fill_discrete(labels=c('segdup better','equal','segdup worse')) +
	theme(legend.title=element_blank()) +
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab("nH") + ylab("Count")
dev.off()

#compare time to Multrec
pdf("figures/nH-sdvmr-time-long.pdf")
plot(means$nH, means$sdTime, col="red", main="Time of segdup vs multrec", xlab="nH", ylim=c(min(means$sdTime-2*sds$sdTime,means$mrTime),max(means$sdTime+2*sds$sdTime,means$mrTime)))
arrows(means$nH, means$sdTime-2*sds$sdTime, means$nH, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3, col="red")
points(means$nH, means$mrTime, col="blue")
arrows(means$nH, means$mrTime-2*sds$mrTime, means$nH, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
legend("topleft", legend=c("segdup","Multrec"), col=c("red","blue"), lty=1)
dev.off()

ggplot(data=means) + geom_point(aes(nH,sdTime,col="red")) + geom_point(aes(nH,mrTime,col="blue")) +
	geom_errorbar(aes(nH,ymin=sdTime-2*sds$sdTime,ymax=sdTime+2*sds$sdTime,col="red")) + 
	geom_errorbar(aes(nH,ymin=mrTime-2*sds$mrTime,ymax=mrTime+2*sds$mrTime,col="blue")) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	scale_colour_manual(values = c('red'='red','blue'='blue'), labels = c('segdup','MultRec'))
	theme(legend.title=element_blank(), legend.position = c(0.05,0.95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
	xlab("nH") + ylab("Time (s)")


#compare cost to true reconciliation
pdf("figures/nH-sdvtr-cost-long.pdf")
plot(means$nH, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="nH", ylim=c(min(means$propTrueCost-2*sds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
arrows(means$nH, means$propTrueCost-2*sds$propTrueCost, means$nH, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
dev.off()


#varying nP
nPdat <- dat[dat$nH==nH & dat$rB==rB & dat$pJ == pJ,]

totals <- aggregate(nPdat, by = list(nPdat$nP), FUN = sum)
means <- aggregate(nPdat, by = list(nPdat$nP), FUN = mean)
sds <- aggregate(nPdat, by = list(nPdat$nP), FUN = sd)/10

#compare cost to Multrec - proportional increase
pdf("figures/nP-sdvmr-cost-long.pdf")
plot(means$nP, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="nP", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
arrows(means$nP, means$propSdCost-2*sds$propSdCost, means$nP, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
dev.off()

#compare cost to Multrec - counts
pdf("figures/nP-sdvmr-counts-long.pdf")
barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$nP, beside=T, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="nP")
legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
dev.off()

#compare time to Multrec
pdf("figures/nP-sdvmr-time-long.pdf")
plot(means$nP, means$sdTime, col="red", main="Time of segdup vs multrec", xlab="nP", ylim=c(min(means$sdTime-2*sds$sdTime,means$mrTime),max(means$sdTime+2*sds$sdTime,means$mrTime)))
arrows(means$nP, means$sdTime-2*sds$sdTime, means$nP, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3, col="red")
points(means$nP, means$mrTime, col="blue")
arrows(means$nP, means$mrTime-2*sds$mrTime, means$nP, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
legend("topleft", legend=c("segdup","Multrec"), col=c("red","blue"), lty=1)
dev.off()

#compare cost to true reconciliation
pdf("figures/nP-sdvtr-cost-long.pdf")
plot(means$nP, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="nP", ylim=c(min(means$propTrueCost-2*sds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
arrows(means$nP, means$propTrueCost-2*sds$propTrueCost, means$nP, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
dev.off()


#varying rB
rBdat <- dat[dat$nH==nH & dat$nP==nP & dat$pJ == pJ,]

totals <- aggregate(rBdat, by = list(rBdat$rB), FUN = sum)
means <- aggregate(rBdat, by = list(rBdat$rB), FUN = mean)
sds <- aggregate(rBdat, by = list(rBdat$rB), FUN = sd)/10

#compare cost to Multrec - proportional increase
pdf("figures/rB-sdvmr-cost-long.pdf")
plot(means$rB, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="rB", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
arrows(means$rB, means$propSdCost-2*sds$propSdCost, means$rB, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
dev.off()

#compare cost to Multrec - counts
pdf("figures/rB-sdvmr-counts-long.pdf")
barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$rB, beside=T, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="rB")
legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
dev.off()

#compare time to Multrec
pdf("figures/rB-sdvmr-time-long.pdf")
plot(means$rB, means$sdTime, col="red", main="Time of segdup vs multrec", xlab="rB", ylim=c(min(means$sdTime-2*sds$sdTime,means$mrTime),max(means$sdTime+2*sds$sdTime,means$mrTime)))
arrows(means$rB, means$sdTime-2*sds$sdTime, means$rB, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3, col="red")
points(means$rB, means$mrTime, col="blue")
arrows(means$rB, means$mrTime-2*sds$mrTime, means$rB, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
legend("topleft", legend=c("segdup","Multrec"), col=c("red","blue"), lty=1)
dev.off()

#compare cost to true reconciliation
pdf("figures/rB-sdvtr-cost-long.pdf")
plot(means$rB, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="rB", ylim=c(min(means$propTrueCost-2*sds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
arrows(means$rB, means$propTrueCost-2*sds$propTrueCost, means$rB, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
dev.off()


#varying pJ
pJdat <- dat[dat$nH == nH & dat$nP==nP & dat$rB==rB,]

totals <- aggregate(pJdat, by = list(pJdat$pJ), FUN = sum)
means <- aggregate(pJdat, by = list(pJdat$pJ), FUN = mean)
sds <- aggregate(pJdat, by = list(pJdat$pJ), FUN = sd)/10

#compare cost to Multrec - proportional increase
pdf("figures/pJ-sdvmr-cost-long.pdf")
plot(means$pJ, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="pJ", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
arrows(means$pJ, means$propSdCost-2*sds$propSdCost, means$pJ, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
dev.off()

#compare cost to Multrec - counts
pdf("figures/pJ-sdvmr-counts-long.pdf")
barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$pJ, beside=T, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="pJ")
legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
dev.off()

#compare time to Multrec
pdf("figures/pJ-sdvmr-time-long.pdf")
plot(means$pJ, means$sdTime, col="red", main="Time of segdup vs multrec", xlab="pJ", ylim=c(min(means$sdTime-2*sds$sdTime,means$mrTime),max(means$sdTime+2*sds$sdTime,means$mrTime)))
arrows(means$pJ, means$sdTime-2*sds$sdTime, means$pJ, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3, col="red")
points(means$pJ, means$mrTime, col="blue")
arrows(means$pJ, means$mrTime-2*sds$mrTime, means$pJ, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
legend("topleft", legend=c("segdup","Multrec"), col=c("red","blue"), lty=1)
dev.off()

#compare cost to true reconciliation
pdf("figures/pJ-sdvtr-cost-long.pdf")
plot(means$pJ, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="pJ", ylim=c(min(means$propTrueCost-2*sds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
arrows(means$pJ, means$propTrueCost-2*sds$propTrueCost, means$pJ, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
dev.off()
