x = read.table("../../../uniform-landscape/gaussian-dispersal/analysis/timings-gaia.txt", header=FALSE)
y = read.table("../../../uniform-landscape/gaussian-dispersal/analysis/timings-wohns.txt", header=FALSE)

x2 = read.csv("../../../uniform-landscape/gaussian-dispersal/analysis/ancestor-estimates.csv", header=FALSE)
y2 = read.csv("../../../uniform-landscape/gaussian-dispersal/ancestor-estimates-wohns.csv", header=FALSE)

e1 = tapply(x2[,5], list(x2[,1],x2[,2]), median)
e2 = tapply(y2[,5], list(y2[,1],y2[,2]), median)

x = c(x[,1], read.table("../timings-gaia.txt", header=FALSE)[,1])
y = c(y[,1], read.table("../timings-wohns.txt", header=FALSE)[,1])

x1 = read.csv("../ancestor-estimates.csv", header=FALSE)
x2 = read.csv("../ancestor-estimates-pareto.csv", header=FALSE)

e1 = c(c(e1), c(tapply(x1[,4], x1[,1], median), tapply(x2[,4], x2[,1], median)))

y1 = read.csv("../ancestor-estimates-wohns.csv", header=FALSE)

idx = 1:802996
e2 = c(c(e2), c(tapply(y1[idx,4], y1[idx,1], median), tapply(y1[-idx,4], y1[-idx,1], median)))


png(file="gaia-wohns-compare.png", res=600, width=9, height=4, units="in")
par(mfrow=c(1,2))
plot(y/x, type='l', ylim=c(0.9,2.5), las=1,
    xlab="Simulation index", ylab="")
points(y/x, cex=0.3,pch=19)
title("x-fold speedup of GAIA compared to Wohns", cex.main=0.8)
abline(h=1,lty=2)
plot(e2/e1, type='l', xlab="Simulation index", las=1,
    ylab="")
points(e2/e1, cex=0.3,pch=19)
title("x-fold accuracy improvement of GAIA compared to Wohns", cex.main=0.8)
abline(h=1,lty=2)
dev.off()