
x = read.csv('../flux/ooa-flux.csv')

y = split(x, x[,1])

par(mfrow=c(1, 3))

plot(0,0,xlim=c(1,n),ylim=c(0,1),type='n',las=1,xlab="index",ylab="flux")
lines(1:n, y$north$flux_north)
points(1:n, y$north$flux_north, pch=19)
lines(1:n, y$north$flux_south,, col=8)
points(1:n, y$north$flux_south, pch=19, col=8)
title("OOA North")
legend("topleft", col=c(1,8), legend=c("flux north", "flux south"), bty='n', lty=c(1,1), cex=0.8)


plot(0,0,xlim=c(1,n),ylim=c(0,1),type='n',las=1,xlab="index",ylab="flux")
lines(1:n, y$south$flux_north)
points(1:n, y$south$flux_north, pch=19)
lines(1:n, y$south$flux_south, col=8)
points(1:n, y$south$flux_south, pch=19, col=8)
title("OOA South")
legend("topleft", col=c(1,8), legend=c("flux north", "flux south"), bty='n', lty=c(1,1), cex=0.8)

plot(0,0,xlim=c(1,n),ylim=c(0,1),type='n',las=1,xlab="index",ylab="flux")
lines(1:n, y$both$flux_north)
points(1:n, y$both$flux_north, pch=19)
lines(1:n, y$both$flux_south, col=8)
points(1:n, y$both$flux_south, pch=19, col=8)
title("OOA North and South")
legend("topleft", col=c(1,8), legend=c("flux north", "flux south"), bty='n', lty=c(1,1), cex=0.8)

dev.print(pdf, file="ooa-flux.pdf")

