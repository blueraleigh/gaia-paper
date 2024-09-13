library(gaia)
library(sf)
library(ape)

mpr = readRDS("../../analysis/results/mpr-chr18p.rds")

ts = treeseq_load("../../data/trees/chr18p.trees")
landgrid = st_read("../../data/landgrid.gpkg")
land = st_read("../../data/land.gpkg")
coords = read.csv("../../data/trees/sample-coords.csv")
sites = st_coordinates(st_centroid(landgrid))

pal = LSD::colorpalette("spectral", rev=FALSE)


treeseq_sample(ts, 1)

phy = treeseq_to_phylo(ts)[[224]]
ancestors = nodepath(phy)
nodes = treeseq_nodes(ts)

tip.id = 186L

bbox = st_bbox(landgrid)
par(xpd=NA)
pmat = persp(
    seq(bbox[1],bbox[3],,10),
    seq(bbox[2],bbox[4],,10),
    matrix(NA,10,10),
    zlim=c(0, 10e6),
    box=FALSE,
    phi=25,
    theta=15,
    scale=TRUE,
    r = sqrt(3), 
    d = 1,
    expand=0.5
)
z = rev(c(12e6,8e6,4e6,0,-4e6))


k=1
for (i in 1:length(landgrid[[1]])) {
    for (ii in 1:length(landgrid[[1]][[i]])) {
        xy = sf:::p_bind(landgrid[[1]][[i]][[ii]])
        xy = trans3d(xy[,1],xy[,2], z[k], pmat)
        polypath(xy$x, xy$y, lwd=0.25, border=8)
    }
}
for (i in 1:length(land[[1]][[1]])) {
    xy = sf:::p_bind(land[[1]][[1]][[i]])
    xy = trans3d(xy[,1],xy[,2], z[k], pmat)
    polypath(xy$x, xy$y, lwd=0.5, border=1)
}
xy = trans3d(
    coords[tail(phy$node.id[ancestors[[tip.id]]]+1,1), 2],
    coords[tail(phy$node.id[ancestors[[tip.id]]]+1,1), 3], z[k], pmat)
points(xy$x, xy$y, pch=21, bg="white", lwd=0.5)
tl = trans3d(bbox[3], bbox[4], z[k], pmat)
text(tl$x, tl$y, "contemporary sample", 
    cex=0.6, pos=3, xpd=NA)

k = 2
for (j in rev(c(1,4,8,14))) {
    a = phy$node.id[ ancestors[[tip.id]][j] ] + 1
    s = 1 - (max(mpr$mpr[,a]) - mpr$mpr[,a]) / (max(mpr$mpr[,a]) - min(mpr$mpr[,a]))
    colv = scales::col_numeric(pal, range(s))(s)
    for (i in 1:length(landgrid[[1]])) {
        for (ii in 1:length(landgrid[[1]][[i]])) {
            xy = sf:::p_bind(landgrid[[1]][[i]][[ii]])
            xy = trans3d(xy[,1],xy[,2], z[k], pmat)
            polypath(xy$x, xy$y, lwd=0.25, col=scales::alpha(colv[i], 1), border=8)
        }
    }
    for (i in 1:length(land[[1]][[1]])) {
        xy = sf:::p_bind(land[[1]][[1]][[i]])
        xy = trans3d(xy[,1],xy[,2], z[k], pmat)
        polypath(xy$x, xy$y, lwd=0.5, border=1)
    }
    xy = trans3d(sites[which.min(s),1], sites[which.min(s),2], z[k], pmat)
    points(xy$x, xy$y, pch=21, bg="white", lwd=0.5)
    tl = trans3d(bbox[3], bbox[4], z[k], pmat)
    if (k == 2)
        text(tl$x, tl$y, sprintf("%.1f kya ancestor", nodes$time[a]*25/1000), 
            cex=0.6, pos=3, xpd=NA)
    else
        text(tl$x, tl$y, sprintf("%.0f kya ancestor", nodes$time[a]*25/1000), 
            cex=0.6, pos=3, xpd=NA)
    k = k+1
}
dev.print(pdf, file="uncertainty.pdf")
