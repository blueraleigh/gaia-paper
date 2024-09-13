cost.mat = data.matrix(read.csv("../../../../hgdp/data/landgrid-distmat.csv", row.names=1))

y = c()
for (i in 1:25) {
    x = read.csv(sprintf("../mpr/ooa-both-mpr-%d.csv", i))
    y = rbind(y, cbind(x[,1], cost.mat[cbind(x[,2],x[,3])] / 22))
}


for (i in 1:25) {
    x = read.csv(sprintf("../mpr/ooa-north-mpr-%d.csv", i))
    y = rbind(y, cbind(x[,1], cost.mat[cbind(x[,2],x[,3])] / 22))
}

for (i in 1:25) {
    x = read.csv(sprintf("../mpr2/ooa-south-mpr-%d.csv", i))
    y = rbind(y, cbind(x[,1], cost.mat[cbind(x[,2],x[,3])] / 22))
}

LSD::heatscatter(y[,1],y[,2], log='x')

f = c()
for (i in 1:25) {
    f = rbind(f, read.csv(sprintf("../flux/ooa-both-flux-%d.csv", i), row.names=1))
}
for (i in 1:25) {
    f = rbind(f, read.csv(sprintf("../flux/ooa-north-flux-%d.csv", i), row.names=1))
}
for (i in 1:25) {
    f = rbind(f, read.csv(sprintf("../flux/ooa-south-flux-%d.csv", i), row.names=1))
}

png(file="ooa-flux-cor.png", width=10, height=5,units="in",res=300)

par(mfrow=c(1,2))

LSD::heatscatter(y[,1],y[,2], log='x', bty='l', main='', las=1,
    xlab='Ancestor age (years)', ylab='Ancestor location error')

LSD::heatscatter(f[,1],f[,2], log='xy', bty='l', xlab='True flux', 
    ylab='Estimated flux', main='', las=1)
abline(0,1)
legend('topleft', bty='n', legend=sprintf('r = %.2f', cor(f[,1],f[,2])))
dev.off()
#dev.print(pdf, file='ooa-flux-cor.pdf')
