library(gaia)

filenames = c()
for (s in seq(0.2,2.0,0.2)) {
    filenames = c(filenames, sprintf("tree-S%.1f-R%d.trees", s, 1:10))
}
par(mfrow=c(3,2), mar=c(2,2,2,2))
j = sample(1:100, 3)
for (i in 1:3)
{

ts = treeseq_load(sprintf('../trees/%s',filenames[j[i]]))

S = as.numeric(substr(strsplit(filenames[j[i]], '-')[[1]][2],2,4))

nodes = treeseq_nodes(ts)

ind = treeseq_individuals(ts)

xy = t(sapply(ind[,2], '[', 1:2))

plot(xy[nodes$individual_id[nodes$time==9999]+1, ], cex=0.3,xlim=c(0,1000*S),
    ylim=c(0,1000*S),xaxt='n', yaxt='n', xlab='', ylab='', main="initial distribution")
plot(xy[nodes$individual_id[nodes$time==0]+1, ], cex=0.05,xlim=c(0,1000*S),
    ylim=c(0,1000*S),xaxt='n', yaxt='n', xlab='', ylab='', main="final distribution")

}

dev.print(pdf, file="homog-init.pdf")
