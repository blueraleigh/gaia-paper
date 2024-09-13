library(gaia)

par(mfrow=c(3,3), mar=c(2,2,2,2))
j = sample(1:20, 3)
for (i in 1:3)
{

ts = treeseq_load(sprintf("../trees/tree-%d.trees", j[i]))
u = data.matrix(read.csv(sprintf("../landscapes/landscape-%d.csv", j[i]),header=FALSE))

nodes = treeseq_nodes(ts)
dim(nodes)

ind = treeseq_individuals(ts)

xy = t(sapply(ind[,2], '[', 1:2))

image(t(u[100:1,]), col=LSD::colorpalette("spectral", rev=TRUE),
    xaxt='n', yaxt='n', main="spatial carrying capacity")
plot(xy[nodes$individual_id[nodes$time==9999]+1, ], cex=0.3, xlim=c(0,25),
    ylim=c(0,25), xaxt='n', yaxt='n', xlab='', ylab='', main="initial distribution")
plot(xy[nodes$individual_id[nodes$time==0]+1, ], cex=0.05, xlim=c(0,25),
    ylim=c(0,25), xaxt='n', yaxt='n', xlab='', ylab='', main="final distribution")

}

dev.print(pdf, file="hetero-init.pdf")