library(ggplot2)
theme_set(theme_bw(base_size=20))

gdat <- read.csv("guigo_temps.csv",header=T)

plot(s7~t,data=gdat)
points(s10~t,data=gdat,col="red")

pdf("guigo_temps.pdf")
with(gdat,plot(s7/(s7+s10) ~ t, xlab="Final temperature", ylab="Proportion of s7 duplications"))
dev.off()

pdf("guigo_temps.pdf")
ggplot(data=gdat) + geom_point(aes(t,s7/(s7+s10))) +
	scale_x_continuous(breaks=seq(0,10,by=2)) +
	xlab("Final temperature") + ylab("Proportion of reconciliations")
dev.off()
