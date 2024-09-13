
x = read.csv("../rate-estimates.csv", header=FALSE)

png(file="rate-estimates.png", width=7, height=7, units="in", res=300)

plot(
x[,4], x[,3],las=1,
ylab=expression(sigma[italic(P)]),xlab=expression(sigma[italic(e)]),
bty="l",type="n",ylim=c(0.1,0.7))
abline(0,1)
points(
x[,4], x[,3], bg="grey",pch=21)

dev.off()

