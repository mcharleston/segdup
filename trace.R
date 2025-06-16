trace <- read.csv("segdup-trace.csv")
par(mfrow=c(1,2))
plot(l ~ i, data = trace, type="l")
plot(d ~ i, data = trace, type="l")
#plot(t ~ i, data = trace, type="l")

library(ggplot2)
theme_set(theme_bw(base_size=20))

pdf("experiments/guigo-trace.pdf")
ggplot(data=trace) + geom_line(aes(i,s)) +
	xlab("Iteration") + ylab("Cost")
dev.off()

pdf("experiments/guigo-dups.pdf")
ggplot(data=trace) + geom_line(aes(i,d)) +
	xlab("Iteration") + ylab("Duplications")
dev.off()
