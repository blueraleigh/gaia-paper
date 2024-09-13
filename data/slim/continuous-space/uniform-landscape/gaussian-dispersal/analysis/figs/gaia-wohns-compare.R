x = read.table("../timings-gaia.txt", header=FALSE)
y = read.table("../timings-wohns.txt", header=FALSE)

x2 = read.csv("../ancestor-estimates-v2.csv", header=FALSE)
y2 = read.csv("../ancestor-estimates-wohns.csv", header=FALSE)

e1 = tapply(x2[,5], list(x2[,1],x2[,2]), median)
e2 = tapply(y2[,5], list(y2[,1],y2[,2]), median)


png(file="gaia-wohns-compare.png", res=600, width=9, height=4, units="in")
par(mfrow=c(1,2))
plot(y[,1]/x[,1], type='l', ylim=c(0.9,2.5), las=1,
    xlab="Simulation index", ylab="")
title("x-fold speedup of GAIA compared to Wohns", cex.main=0.8)
abline(h=1,lty=2)
plot(c(e2/e1), type='l', xlab="Simulation index", las=1,
    ylab="")
title("x-fold accuracy improvement of GAIA compared to Wohns", cex.main=0.8)
abline(h=1,lty=2)
dev.off()