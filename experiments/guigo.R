gdat <- read.csv("guigo_temps.csv",header=T)

plot(s7~t,data=gdat)
points(s10~t,data=gdat,col="red")

pdf("guigo_temps.pdf")
with(gdat,plot(s7/(s7+s10) ~ t, xlab="Final temperature", ylab="Proportion of s7 duplications"))
dev.off()
