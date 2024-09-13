library(sf)

data = read.csv("../../../../hgdp/data/trees/sample-states.csv")
pops_to_sample = sort(unique(data[,2]))

landgrid = st_read("../../../../hgdp/data/landgrid.gpkg")

y = c()
for (i in 1:25) {
    x = read.csv(sprintf("../mpr/ooa-both-mpr-%d.csv", i))
    idx = which(x$node_time == 19999)
    y = rbind(y, x[idx,2:3])
}


for (i in 1:25) {
    x = read.csv(sprintf("../mpr/ooa-north-mpr-%d.csv", i))
    idx = which(x$node_time == 19999)
    y = rbind(y, x[idx,2:3])
}


for (i in 1:25) {
    x = read.csv(sprintf("../mpr/ooa-south-mpr-%d.csv", i))
    idx = which(x$node_time == 19999)
    y = rbind(y, x[idx,2:3])
}


par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(landgrid$geom, border=8, lwd=0.5)
points(st_coordinates(st_centroid(landgrid$geom))[pops_to_sample,], 
    cex=0.8, lwd=0.5)
points(st_coordinates(st_centroid(landgrid$geom))[58,,drop=FALSE], pch=20, 
    col=2, cex=0.5)

legend("top", legend=c("modern day samples", "founding individuals"),
    bty="n", pch=c(1, 20), col=c(1,2), ncol=2, cex=0.8, inset=.2, pt.lwd=0.5)


plot(landgrid$geom, border=8, lwd=0.5)

points(st_coordinates(st_centroid(landgrid$geom)), lwd=0.5,
    cex=0.15*sqrt(tabulate(y[,2], 177)), pch=21, bg='white')

points(st_coordinates(st_centroid(landgrid$geom))[58,,drop=FALSE], pch=20, 
    col=2, cex=0.5)

legend("top", legend="GAIA estimates of founding individuals",
    bty="n", pch=1, col=1, ncol=1, pt.cex=1, cex=0.8, inset=.2, pt.lwd=0.5)


dev.print(pdf, 'ooa-ancestor-estimates.pdf')
