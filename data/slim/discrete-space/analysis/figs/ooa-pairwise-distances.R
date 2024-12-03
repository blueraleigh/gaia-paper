adjmat = read.csv("../../simulations/landgrid-adjmat-ooa-both.csv",header=TRUE)
A = data.matrix(adjmat)
dimnames(A) = NULL
D = igraph::distances(igraph::graph_from_adjacency_matrix(A,mode="undirected"))

par(mfrow=c(1,4))
ages = c(500, 1000, 5000, 10000, 20000)
for(k in 1:4)
{
d = matrix(0, 23, 23)
for (i in 1:20)
{
x = read.csv(sprintf("../mpr/ooa-both-mpr-%d.csv",i))
idx = which(x[,1] < ages[k+1] & x[,1] >= ages[k])
# 1. randomly sample 100 nodes from the tree sequence
#    1000 times (w/o replacement)
# 2. for each sample, compute the estimated distance between every
#    pair of nodes in the sample
# 3. for each sample, create a frequency table mapping true pairwise
#    distances to estimated pairwise distances
# 4. average these tables over all samples
# 5. repeat for each simulation replicate
for (j in 1:1000)
{
s = sample(idx, 100)
g = data.matrix(expand.grid(s, s))
g = g[g[,1] != g[,2],]

true_path_distance = D[cbind(x[g[,1],2], x[g[,2],2])]
estimated_path_distance = D[cbind(x[g[,1],3], x[g[,2],3])]

true_path_distance = factor(true_path_distance, levels=0:22)
estimated_path_distance = factor(estimated_path_distance, levels=0:22)

d = d + table(true_path_distance, estimated_path_distance) / 20000
}
}
image(0:22, 0:22, d, col=LSD::colorpalette("heat", nrcol=264),las=1,
    xlab="True pairwise distance", ylab="Estimated pairwise distance")
abline(0,1)
title(sprintf("%d - %d years ago", ages[k], ages[k+1]), cex=0.8, adj=0)
text(6, 21.5, "Avg number of pairs", cex=0.8)
rect(seq(1, 12, , 265)[-265], 20, seq(1, 12, , 265)[-1], 21,
    col=LSD::colorpalette("heat", nrcol=264), border=NA)
text(seq(1, 12, , 265)[2], 19.5, 1, cex=0.6)
text(seq(1, 12, , 265)[264], 19.5, round(max(d)), cex=0.6)
}

#i = 12
#colv = rep(8,23)
#colv[i] = 1
#barplot(d[i,], col=colv, las=1)
#mtext("Estimated pairwise distance", 1, line=3)
#text(i+1.5, max(d[i,]), "True pairwise distance", pos=4, xpd=NA)
