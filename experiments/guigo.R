gdat <- read.csv("guigo_temps.csv",header=T)

plot(s7~t,data=gdat)
points(s10~t,data=gdat,col="red")

with(gdat,plot(s7/(s7+s10)))
