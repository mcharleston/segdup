library(ggplot2)
library(reshape)
theme_set(theme_bw(base_size=20))

#system("python3 process.py")
dat <- read.csv("allresults.csv", header=T)
datlong <- read.csv("allresults-long.csv", header=T)

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

datlong$trueCost <- d*datlong$nAllDupEvents+l*datlong$nLineageSort

datlong$propSdCost <- datlong$cost/datlong$mrCost
datlong$propTrueCost <- datlong$cost/datlong$trueCost

datlong$sdBetter <- datlong$cost < datlong$mrCost
datlong$same <- datlong$cost == datlong$mrCost
datlong$mrBetter <- datlong$mrCost < datlong$cost


#varying nH
nHdat <- dat[dat$nP==nP & dat$rB==rB & dat$pJ == pJ,]
nHdatlong <- datlong[datlong$nP==nP & datlong$rB==rB & datlong$pJ == pJ,]

totals <- aggregate(nHdat, by = list(nHdat$nH), FUN = sum)
means <- aggregate(nHdat, by = list(nHdat$nH), FUN = mean)
sds <- aggregate(nHdat, by = list(nHdat$nH), FUN = sd)/10

ltotals <- aggregate(nHdatlong, by = list(nHdatlong$nH), FUN = sum)
lmeans <- aggregate(nHdatlong, by = list(nHdatlong$nH), FUN = mean)
lsds <- aggregate(nHdatlong, by = list(nHdatlong$nH), FUN = sd)/10

#compare cost to Multrec - proportional increase
#pdf("figures/nH-sdvmr-cost-long.pdf")
#plot(means$nH, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="nH", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
#arrows(means$nH, means$propSdCost-2*sds$propSdCost, means$nH, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
#dev.off()

pdf("figures/nH-sdvmr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(nH-1,propSdCost,col="black")) +
	geom_errorbar(aes(nH-1,ymin=propSdCost-2*sds$propSdCost,ymax=propSdCost+2*sds$propSdCost),width=2) + 
	geom_point(aes(nH,lmeans$propSdCost,col="red")) +
	geom_errorbar(aes(nH,ymin=lmeans$propSdCost-2*lsds$propSdCost,ymax=lmeans$propSdCost+2*lsds$propSdCost),col="red",width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(n[H])) + ylab("Proportional cost")
dev.off()

#compare cost to Multrec - counts
#pdf("figures/nH-sdvmr-counts-long.pdf")
#barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$nH, beside=T, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="nH")
#legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
#dev.off()

pdf("figures/nH-sdvmr-counts-both.pdf")
tm <- melt(data.frame(nH=means$nH,sd=means$sdBetter*100,mr=means$mrBetter*100,lsd=lmeans$sdBetter*100,lmr=lmeans$mrBetter*100), id.vars=1)
ggplot(data=tm) + geom_bar(aes(nH, value, fill=variable), position="dodge", stat="identity") +
	scale_fill_manual(values=alpha(c('red','blue','red','blue'),c(0.25,0.25,1,1)), labels=c(bquote(segdup~better~(10^4)),bquote(segdup~worse~(10^4)),bquote(segdup~better~(10^5)),bquote(segdup~worse~(10^5)))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	xlab(bquote(n[H])) + ylab("Count")
dev.off()

#compare time to Multrec
#pdf("figures/nH-sdvmr-time-long.pdf")
#plot(means$nH, means$sdTime, col="red", main="Time of segdup vs multrec", xlab="nH", ylim=c(min(means$sdTime-2*sds$sdTime,means$mrTime),max(means$sdTime+2*sds$sdTime,means$mrTime)))
#arrows(means$nH, means$sdTime-2*sds$sdTime, means$nH, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3, col="red")
#points(means$nH, means$mrTime, col="blue")
#arrows(means$nH, means$mrTime-2*sds$mrTime, means$nH, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
#legend("topleft", legend=c("segdup","Multrec"), col=c("red","blue"), lty=1)
#dev.off()

pdf("figures/nH-sdvmr-time-both.pdf")
ggplot(data=means) + geom_point(aes(nH-1,sdTime,col="black")) + 
	geom_errorbar(aes(nH-1,ymin=sdTime-2*sds$sdTime,ymax=sdTime+2*sds$sdTime,col="black"),width=2) + 
	geom_point(aes(nH,lmeans$sdTime,col="red")) + 
	geom_errorbar(aes(nH,ymin=lmeans$sdTime-2*lsds$sdTime,ymax=lmeans$sdTime+2*lsds$sdTime,col="red"),width=2) +
	geom_point(aes(nH+1,(mrTime+lmeans$mrTime)/2,colour="blue")) +
	geom_errorbar(aes(nH+1,ymin=pmax(0,(mrTime+lmeans$mrTime)/2-2*(sds$mrTime+lsds$mrTime)/2/sqrt(2)),ymax=(mrTime+lmeans$mrTime)/2+2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),colour="blue"),width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red','blue'), breaks=c('black','red','blue'),labels = c(bquote(segdup~(10^4)),bquote(segdup~(10^5)),'MultRec')) +
	theme(legend.title=element_blank(), legend.position = c(0.95,0.05), legend.justification = c("right", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
	#coord_cartesian(ylim=c(0,max(lmeans$sdTime+2*lsds$sdTime))) +
	scale_y_continuous(trans='log',breaks=10^(-1:3),labels=c("0.1","1","10","100","1000")) +
	xlab(bquote(n[H])) + ylab("Time (s)")
dev.off()

#compare cost to true reconciliation
#pdf("figures/nH-sdvtr-cost-long.pdf")
#plot(means$nH, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="nH", ylim=c(min(means$propTrueCost-2*sds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
#arrows(means$nH, means$propTrueCost-2*sds$propTrueCost, means$nH, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
#dev.off()

pdf("figures/nH-sdvtr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(nH-1,propTrueCost,col="black")) +
	geom_errorbar(aes(nH-1,ymin=propTrueCost-2*sds$propTrueCost,ymax=propTrueCost+2*sds$propTrueCost,col="black"),width=2) +
	geom_point(aes(nH,lmeans$propTrueCost,col="red")) +
	geom_errorbar(aes(nH,ymin=lmeans$propTrueCost-2*lsds$propTrueCost,ymax=lmeans$propTrueCost+2*lsds$propTrueCost,col="red"),width=2) +
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .05), legend.justification = c("right", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(n[H])) + ylab("Proportional cost")
dev.off()

#varying nP
nPdat <- dat[dat$nH==nH & dat$rB==rB & dat$pJ == pJ,]
nPdatlong <- datlong[datlong$nH==nH & datlong$rB==rB & datlong$pJ == pJ,]

totals <- aggregate(nPdat, by = list(nPdat$nP), FUN = sum)
means <- aggregate(nPdat, by = list(nPdat$nP), FUN = mean)
sds <- aggregate(nPdat, by = list(nPdat$nP), FUN = sd)/10

ltotals <- aggregate(nPdatlong, by = list(nPdatlong$nP), FUN = sum)
lmeans <- aggregate(nPdatlong, by = list(nPdatlong$nP), FUN = mean)
lsds <- aggregate(nPdatlong, by = list(nPdatlong$nP), FUN = sd)/10

#compare cost to Multrec - proportional increase
#pdf("figures/nP-sdvmr-cost-long.pdf")
#plot(means$nP, means$propSdCost, main="Proportional increase in cost of segdup vs multrec", xlab="nP", ylim=c(min(means$propSdCost-2*sds$propSdCost),max(means$propSdCost+2*sds$propSdCost)))
#arrows(means$nP, means$propSdCost-2*sds$propSdCost, means$nP, means$propSdCost+2*sds$propSdCost, length=0.05, angle=90, code=3)
#dev.off()

pdf("figures/nP-sdvmr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(nP-1,propSdCost,col="black")) +
	geom_errorbar(aes(nP-1,ymin=propSdCost-2*sds$propSdCost,ymax=propSdCost+2*sds$propSdCost),width=2) + 
	geom_point(aes(nP,lmeans$propSdCost,col="red")) +
	geom_errorbar(aes(nP,ymin=lmeans$propSdCost-2*lsds$propSdCost,ymax=lmeans$propSdCost+2*lsds$propSdCost),col="red",width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(n[P])) + ylab("Proportional cost")
dev.off()

#compare cost to Multrec - counts
#pdf("figures/nP-sdvmr-counts-long.pdf")
#
#barplot(t(as.matrix(totals[,c("sdBetter","same","mrBetter")])), names.arg=means$nP, beside=F, col=c("red","green","blue"), main="Counts where segdup is better/equal/worse than multrec", xlab="nP")
#legend("topright", legend=c("segdup better","equal","segdup worse"), fill=c("red","green","blue"))
#
#dev.off()
#
#pdf("figures/nP-sdvmr-counts-both.pdf")
#barplot(t(as.matrix(cbind(totals[,c("sdBetter","mrBetter")],ltotals[,c("sdBetter","mrBetter")]))), names.arg=means$nP, beside=T, col=c("red","blue"), density = c(50,50,1000,1000), main="Counts where segdup is better/equal/worse than multrec", xlab="nP")
#legend("topleft", legend=c("segdup better (10^4)","segdup worse (10^4)","segdup better (10^5)","segdup worse (10^5)"), fill=c("red","blue"), density= c(50,50,1000,1000))
#dev.off()

pdf("figures/nP-sdvmr-counts-both.pdf")
tm <- melt(data.frame(nP=means$nP,sd=means$sdBetter*100,mr=means$mrBetter*100,lsd=lmeans$sdBetter*100,lmr=lmeans$mrBetter*100), id.vars=1)
ggplot(data=tm) + geom_bar(aes(nP, value, fill=variable), position="dodge", stat="identity") +
	scale_fill_manual(values=alpha(c('red','blue','red','blue'),c(0.25,0.25,1,1)), labels=c(bquote(segdup~better~(10^4)),bquote(segdup~worse~(10^4)),bquote(segdup~better~(10^5)),bquote(segdup~worse~(10^5)))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	xlab(bquote(n[P])) + ylab("Count")
dev.off()

#compare time to Multrec
#pdf("figures/nP-sdvmr-time-both.pdf")
#plot(means$nP, means$sdTime, main="Time of segdup vs multrec", xlab="nP", ylim=c(0,max(lmeans$sdTime+2*lsds$sdTime,lmeans$mrTime)))
#arrows(means$nP, means$sdTime-2*sds$sdTime, means$nP, means$sdTime+2*sds$sdTime, length=0.05, angle=90, code=3)
#points(lmeans$nP, lmeans$sdTime, col="red")
#arrows(lmeans$nP, lmeans$sdTime-2*lsds$sdTime, lmeans$nP, lmeans$sdTime+2*lsds$sdTime, length=0.05, angle=90, code=3, col="red")
#points(means$nP, means$mrTime, col="blue")
#arrows(means$nP, means$mrTime-2*sds$mrTime, means$nP, means$mrTime+2*sds$mrTime, length=0.05, angle=90, code=3, col="blue")
#legend("topleft", legend=c("segdup (10^4)","segdup (10^5)","Multrec"), col=c("black","red","blue"), lty=1)
#dev.off()

pdf("figures/nP-sdvmr-time-both.pdf")
ggplot(data=means) + geom_point(aes(nP-1,sdTime,col="black")) + 
	geom_errorbar(aes(nP-1,ymin=sdTime-2*sds$sdTime,ymax=sdTime+2*sds$sdTime,col="black"),width=2) + 
	geom_point(aes(nP,lmeans$sdTime,col="red")) + 
	geom_errorbar(aes(nP,ymin=lmeans$sdTime-2*lsds$sdTime,ymax=lmeans$sdTime+2*lsds$sdTime,col="red"),width=2) +
	geom_point(aes(nP+1,(mrTime+lmeans$mrTime)/2,colour="blue")) +
	geom_errorbar(aes(nP+1,ymin=pmax(0,(mrTime+lmeans$mrTime)/2-2*(sds$mrTime+lsds$mrTime)/2/sqrt(2)),ymax=(mrTime+lmeans$mrTime)/2+2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),colour="blue"),width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red','blue'), breaks=c('black','red','blue'),labels = c(bquote(segdup~(10^4)),bquote(segdup~(10^5)),'MultRec')) +
	theme(legend.title=element_blank(), legend.position = c(0.05,0.95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
	#coord_cartesian(ylim=c(0,1.5*max(lmeans$sdTime+2*lsds$sdTime))) +
	scale_y_continuous(trans='log',breaks=10^(-1:3),labels=c("0.1","1","10","100","1000")) +
	xlab(bquote(n[P])) + ylab("Time (s)")
dev.off()

#compare cost to true reconciliation
#pdf("figures/nP-sdvtr-cost-both.pdf")
#plot(means$nP, means$propTrueCost, main="Proportional decrease in cost of segdup vs true", xlab="nP", ylim=c(min(lmeans$propTrueCost-2*lsds$propTrueCost),max(means$propTrueCost+2*sds$propTrueCost)))
#arrows(means$nP, means$propTrueCost-2*sds$propTrueCost, means$nP, means$propTrueCost+2*sds$propTrueCost, length=0.05, angle=90, code=3)
#points(lmeans$nP, lmeans$propTrueCost,col="red")
#arrows(lmeans$nP, lmeans$propTrueCost-2*lsds$propTrueCost, lmeans$nP, lmeans$propTrueCost+2*lsds$propTrueCost, length=0.05, angle=90, code=3, col="red")
#legend("topright", legend=c("segdup (10^4)","segdup (10^5)"), col=c("black","red"), lty=1)
#dev.off()

pdf("figures/nP-sdvtr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(nP-1,propTrueCost,col="black")) +
	geom_errorbar(aes(nP-1,ymin=propTrueCost-2*sds$propTrueCost,ymax=propTrueCost+2*sds$propTrueCost,col="black"),width=2) +
	geom_point(aes(nP,lmeans$propTrueCost,col="red")) +
	geom_errorbar(aes(nP,ymin=lmeans$propTrueCost-2*lsds$propTrueCost,ymax=lmeans$propTrueCost+2*lsds$propTrueCost,col="red"),width=2) +
	scale_x_continuous(breaks=seq(0,100,by=10)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	coord_cartesian(ylim=c(0.9,0.96)) +
	xlab(bquote(n[P])) + ylab("Proportional cost")
dev.off()



#varying rB
rBdat <- dat[dat$nH==nH & dat$nP==nP & dat$pJ == pJ,]
rBdatlong <- datlong[datlong$nH==nH & datlong$nP==nP & datlong$pJ == pJ,]

totals <- aggregate(rBdat, by = list(rBdat$rB), FUN = sum)
means <- aggregate(rBdat, by = list(rBdat$rB), FUN = mean)
sds <- aggregate(rBdat, by = list(rBdat$rB), FUN = sd)/10

ltotals <- aggregate(rBdatlong, by = list(rBdatlong$rB), FUN = sum)
lmeans <- aggregate(rBdatlong, by = list(rBdatlong$rB), FUN = mean)
lsds <- aggregate(rBdatlong, by = list(rBdatlong$rB), FUN = sd)/10

#compare cost to Multrec - proportional increase
pdf("figures/rB-sdvmr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(rB-0.05,propSdCost,col="black")) +
	geom_errorbar(aes(rB-0.05,ymin=propSdCost-2*sds$propSdCost,ymax=propSdCost+2*sds$propSdCost),width=0.1) + 
	geom_point(aes(rB,lmeans$propSdCost,col="red")) +
	geom_errorbar(aes(rB,ymin=lmeans$propSdCost-2*lsds$propSdCost,ymax=lmeans$propSdCost+2*lsds$propSdCost),col="red",width=0.1) + 
	scale_x_continuous(breaks=seq(1,5,by=1)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(r[B])) + ylab("Proportional cost")
dev.off()

#compare cost to Multrec - counts
pdf("figures/rB-sdvmr-counts-both.pdf")
tm <- melt(data.frame(rB=means$rB,sd=means$sdBetter*100,mr=means$mrBetter*100,lsd=lmeans$sdBetter*100,lmr=lmeans$mrBetter*100), id.vars=1)
ggplot(data=tm) + geom_bar(aes(rB, value, fill=variable), position="dodge", stat="identity") +
	scale_fill_manual(values=alpha(c('red','blue','red','blue'),c(0.25,0.25,1,1)), labels=c(bquote(segdup~better~(10^4)),bquote(segdup~worse~(10^4)),bquote(segdup~better~(10^5)),bquote(segdup~worse~(10^5)))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	scale_x_continuous(breaks=seq(1,5,by=1)) +
	xlab(bquote(r[B])) + ylab("Count")
dev.off()

#compare time to Multrec
pdf("figures/rB-sdvmr-time-both.pdf")
ggplot(data=means) + geom_point(aes(rB-0.05,sdTime,col="black")) + 
	geom_errorbar(aes(rB-0.05,ymin=sdTime-2*sds$sdTime,ymax=sdTime+2*sds$sdTime,col="black"),width=0.1) + 
	geom_point(aes(rB,lmeans$sdTime,col="red")) + 
	geom_errorbar(aes(rB,ymin=lmeans$sdTime-2*lsds$sdTime,ymax=lmeans$sdTime+2*lsds$sdTime,col="red"),width=0.1) +
	geom_point(aes(rB+0.05,(mrTime+lmeans$mrTime)/2,colour="blue")) +
	geom_errorbar(aes(rB+0.05,ymin=(mrTime+lmeans$mrTime)/2-2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),ymax=(mrTime+lmeans$mrTime)/2+2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),colour="blue"),width=0.1) + 
	scale_x_continuous(breaks=seq(1,5,by=1)) +
	scale_colour_manual(values = c('black','red','blue'), breaks=c('black','red','blue'),labels = c(bquote(segdup~(10^4)),bquote(segdup~(10^5)),'MultRec')) +
	theme(legend.title=element_blank(), legend.position = c(0.05,0.95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
	#coord_cartesian(ylim=c(0,max(lmeans$mrTime))) +
	scale_y_continuous(trans='log',breaks=10^(-1:3),labels=c("0.1","1","10","100","1000")) +
	xlab(bquote(r[B])) + ylab("Time (s)")
dev.off()

#compare cost to true reconciliation
pdf("figures/rB-sdvtr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(rB-0.05,propTrueCost,col="black")) +
	geom_errorbar(aes(rB-0.05,ymin=propTrueCost-2*sds$propTrueCost,ymax=propTrueCost+2*sds$propTrueCost,col="black"),width=0.1) +
	geom_point(aes(rB,lmeans$propTrueCost,col="red")) +
	geom_errorbar(aes(rB,ymin=lmeans$propTrueCost-2*lsds$propTrueCost,ymax=lmeans$propTrueCost+2*lsds$propTrueCost,col="red"),width=0.1) +
	scale_x_continuous(breaks=seq(1,5,by=1)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .05), legend.justification = c("right", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(r[B])) + ylab("Proportional cost")
dev.off()


#varying pJ
pJdat <- dat[dat$nH == nH & dat$nP==nP & dat$rB==rB,]
pJdatlong <- datlong[datlong$nH == nH & datlong$nP==nP & datlong$rB==rB,]

totals <- aggregate(pJdat, by = list(pJdat$pJ), FUN = sum)
means <- aggregate(pJdat, by = list(pJdat$pJ), FUN = mean)
sds <- aggregate(pJdat, by = list(pJdat$pJ), FUN = sd)/10

ltotals <- aggregate(pJdatlong, by = list(pJdatlong$pJ), FUN = sum)
lmeans <- aggregate(pJdatlong, by = list(pJdatlong$pJ), FUN = mean)
lsds <- aggregate(pJdatlong, by = list(pJdatlong$pJ), FUN = sd)/10

#compare cost to Multrec - proportional increase
pdf("figures/pJ-sdvmr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(pJ-0.01,propSdCost,col="black")) +
	geom_errorbar(aes(pJ-0.01,ymin=propSdCost-2*sds$propSdCost,ymax=propSdCost+2*sds$propSdCost),width=0.02) + 
	geom_point(aes(pJ,lmeans$propSdCost,col="red")) +
	geom_errorbar(aes(pJ,ymin=lmeans$propSdCost-2*lsds$propSdCost,ymax=lmeans$propSdCost+2*lsds$propSdCost),col="red",width=0.02) + 
	scale_x_continuous(breaks=seq(0,1,by=0.2)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(p[J])) + ylab("Proportional cost")
dev.off()

#compare cost to Multrec - counts
pdf("figures/pJ-sdvmr-counts-both.pdf")
tm <- melt(data.frame(pJ=means$pJ,sd=means$sdBetter*100,mr=means$mrBetter*100,lsd=lmeans$sdBetter*100,lmr=lmeans$mrBetter*100), id.vars=1)
ggplot(data=tm) + geom_bar(aes(pJ, value, fill=variable), position="dodge", stat="identity") +
	scale_fill_manual(values=alpha(c('red','blue','red','blue'),c(0.25,0.25,1,1)), labels=c(bquote(segdup~better~(10^4)),bquote(segdup~worse~(10^4)),bquote(segdup~better~(10^5)),bquote(segdup~worse~(10^5)))) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	scale_x_continuous(breaks=seq(0,1,by=0.2)) +
	xlab(bquote(p[J])) + ylab("Count")
dev.off()

#compare time to Multrec
pdf("figures/pJ-sdvmr-time-both.pdf")
ggplot(data=means) + geom_point(aes(pJ-0.01,sdTime,col="black")) + 
	geom_errorbar(aes(pJ-0.01,ymin=sdTime-2*sds$sdTime,ymax=sdTime+2*sds$sdTime,col="black"),width=0.02) + 
	geom_point(aes(pJ,lmeans$sdTime,col="red")) + 
	geom_errorbar(aes(pJ,ymin=lmeans$sdTime-2*lsds$sdTime,ymax=lmeans$sdTime+2*lsds$sdTime,col="red"),width=0.02) +
	geom_point(aes(pJ+0.01,(mrTime+lmeans$mrTime)/2,colour="blue")) +
	geom_errorbar(aes(pJ+0.01,ymin=(mrTime+lmeans$mrTime)/2-2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),ymax=(mrTime+lmeans$mrTime)/2+2*(sds$mrTime+lsds$mrTime)/2/sqrt(2),colour="blue"),width=0.02) + 
	scale_x_continuous(breaks=seq(0,1,by=0.2)) +
	scale_colour_manual(values = c('black','red','blue'), breaks=c('black','red','blue'),labels = c(bquote(segdup~(10^4)),bquote(segdup~(10^5)),'MultRec')) +
	theme(legend.title=element_blank(), legend.position = c(0.05,0.95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) + 
	#coord_cartesian(ylim=c(0,max(lmeans$mrTime))) +
	scale_y_continuous(trans='log',breaks=10^(-1:3),labels=c("0.1","1","10","100","1000")) +
	xlab(bquote(p[J])) + ylab("Time (s)")
dev.off()

#compare cost to true reconciliation
pdf("figures/pJ-sdvtr-cost-both.pdf")
ggplot(data=means) + geom_point(aes(pJ-0.01,propTrueCost,col="black")) +
	geom_errorbar(aes(pJ-0.01,ymin=propTrueCost-2*sds$propTrueCost,ymax=propTrueCost+2*sds$propTrueCost,col="black"),width=0.02) +
	geom_point(aes(pJ,lmeans$propTrueCost,col="red")) +
	geom_errorbar(aes(pJ,ymin=lmeans$propTrueCost-2*lsds$propTrueCost,ymax=lmeans$propTrueCost+2*lsds$propTrueCost,col="red"),width=0.02) +
	scale_x_continuous(breaks=seq(0,1,by=0.2)) +
	scale_colour_manual(values = c('black','red'), labels = c(bquote(10^4~' iterations'),bquote(10^5~' iterations'))) +
	theme(legend.title=element_blank(), legend.position = c(.95, .05), legend.justification = c("right", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(p[J])) + ylab("Proportional cost")
dev.off()
