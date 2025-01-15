library(sf)


data = read.csv("../../data/trees/sample-states.csv")
landgrid = st_read("../../data/landgrid.gpkg")
land = st_read("../../data/land.gpkg")
site_coords = st_coordinates(st_centroid(landgrid))

flux = readRDS("../../analysis/results/flux-thru-time-subset-avg-chr18p.rds")
Q = readRDS("../../analysis/results/ancestry-thru-time-chr18p.rds")
Z = matrix(0, 177, 201)
i = 1
invisible(apply(Q, 1, function(a) {
    tmp = attr(Q,'scale') * a
    Z[i,] <<- colSums(tmp)
    i <<- i + 1
}, simplify=FALSE))
Z = sweep(Z,2,colSums(Z),'/')

sites = st_centroid(landgrid$geom)

site_coords = st_coordinates(sites)

times = seq(0, 20000, 100)

pdf(file="ooa-flux.pdf", width=7.25, height=4.25, colormodel="cmyk")
l = layout(matrix(c(1,1,1,2,3,4), 3, 2))
#layout.show(l)

par(mar=c(0,0,0,0), oma=c(5,4,2,0), bty="n", xpd=NA, colormodel="cmyk")
plot(colSums(flux[38,69:70,])-colSums(flux[69:70,38,]), times[-201],
    xlim=c(-0.02, 0.02), type='n', las=1, xlab="", ylab="", yaxt="n",
    xaxt="n", bty="n")
segments(0, 0, 0, 20000, col=8)
t = c(200, 50, 25)
se = attr(flux, 'std.error')
for (i in 1:200)
{
    d = sum(flux[38,69:70,i])+flux[67,68,i]
    v = (sum(se[38,69:70,i])+se[67,68,i])*sqrt(100)
    u = d + v
    l = d - v
    segments(max(l, 0), times[i], u, times[i],col="#d8b365")
    points(d, times[i], pch=19, cex=0.2)
    d = sum(flux[69:70,38,i])+flux[68,67,i]
    v = (sum(se[69:70,38,i])+se[68,67,i])*sqrt(100)
    u = d + v
    l = d - v
    segments(-u, times[i], min(-l, 0), times[i],col="#5ab4ac")
    points(-d, times[i], pch=19, cex=0.2)
}
axis(2, las=1, tcl=-0.2, mgp=c(3,0,-0.5), at=c(0,5000,10000,15000,20000),
    labels=c(0,125,250,375,500), cex.axis=0.75)
mtext("Thousands of years ago", 2, line=1.5, cex=0.75)
axis(1, tcl=-0.2, mgp=c(3,0,-0.25), at=c(-.02,-.01,0,.01,.02), labels=c(.02,.01,0,.01,.02), cex.axis=0.75)
mtext("Into Africa", 1, at=-0.02, line=1.5, cex=0.75, adj=0)
mtext("Out of Africa", 1, at=0.02, line=1.5, cex=0.75, adj=1)
mtext("Migration", 1, line=2, cex=0.75)

tmp = Z[,t[1]]
plot(0,0,type='n',xlim=st_bbox(landgrid)[c(2,4)],ylim=st_bbox(landgrid)[c(2,4)],
    xaxt="n", yaxt="n",bty="n",xlab="",ylab="",asp=1)
plot(landgrid$geom,add=TRUE,
    col=scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
        c(0,max(tmp)), na.color="light grey")(tmp), 
    border="light grey", lwd=0.5)
for (j in 1:177)
{
    nbrs = which(flux[,j,t[1]] > 0)
    if (!length(nbrs)) next
    wt = flux[nbrs,j,t[1]]
    nbrs = nbrs[which.max(wt)]
    wt = max(wt)
    m = colSums(sweep(site_coords[nbrs,,drop=FALSE], 1, wt, "*")) / sum(wt)
    arrows(
        0.85*m[1] + 0.15*site_coords[j,1], 
        0.85*m[2] + 0.15*site_coords[j,2], 
        0.85*site_coords[j,1] + 0.15*m[1], 
        0.85*site_coords[j,2] + 0.15*m[2], length=0.025,
        col=scales::alpha(1, 1000*sum(wt)), lwd=0.5)#lwd=sum(wt)*100)
}
title(sprintf("%s kya", times[t[1]]*25/1000), adj=0.1, cex.main=10/12)

tmp = Z[,t[2]]
plot(0,0,type='n',xlim=st_bbox(landgrid)[c(2,4)],ylim=st_bbox(landgrid)[c(2,4)],
    xaxt="n", yaxt="n",bty="n",xlab="",ylab="",asp=1)
plot(landgrid$geom,add=TRUE,
    col=scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
        c(0,max(tmp)), na.color="light grey")(tmp), 
    border="light grey", lwd=0.5)
for (j in 1:177)
{
    nbrs = which(flux[,j,t[2]] > 0)
    if (!length(nbrs)) next
    wt = flux[nbrs,j,t[2]]
    nbrs = nbrs[which.max(wt)]
    wt = max(wt)
    m = colSums(sweep(site_coords[nbrs,,drop=FALSE], 1, wt, "*")) / sum(wt)
    arrows(
        0.85*m[1] + 0.15*site_coords[j,1], 
        0.85*m[2] + 0.15*site_coords[j,2], 
        0.85*site_coords[j,1] + 0.15*m[1], 
        0.85*site_coords[j,2] + 0.15*m[2], length=0.025,
        col=scales::alpha(1, 1000*sum(wt)), lwd=0.5)#lwd=sum(wt)*100)
}
title(sprintf("%s kya", times[t[2]]*25/1000), adj=0.1, cex.main=10/12)
g = st_make_grid(landgrid$geom, square=FALSE, flat_topped=TRUE,
    n=20)

plot(g[rev(c(131,155,179,203,227))], add=TRUE, border="light grey",
    col=scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
        c(0,1), na.color="light grey")(c(0.05, 0.25, 0.45, 0.65, 0.85)),
        lwd=0.5)
text(
    st_coordinates(st_centroid(g[227]))[1],
    st_coordinates(st_centroid(g[227]))[2],
    "Low ancestry", pos=4, cex=0.75)
text(
    st_coordinates(st_centroid(g[131]))[1],
    st_coordinates(st_centroid(g[131]))[2],
    "High ancestry", pos=4, cex=0.75)

plot(g[83], add=TRUE, border="light grey")
arrows(
    st_coordinates(st_centroid(g[95]))[1],
    st_coordinates(st_centroid(g[95]))[2],
    st_coordinates(st_centroid(g[83]))[1],
    st_coordinates(st_centroid(g[83]))[2],
    length=0.025, lwd=0.5
)
text(
    st_coordinates(st_centroid(g[83]))[1],
    st_coordinates(st_centroid(g[83]))[2],
    "Net migration", pos=4, cex=0.75)

tmp = Z[,t[3]]
plot(0,0,type='n',xlim=st_bbox(landgrid)[c(2,4)],ylim=st_bbox(landgrid)[c(2,4)],
    xaxt="n", yaxt="n",bty="n",xlab="",ylab="",asp=1)
plot(landgrid$geom,add=TRUE,
    col=scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
        c(0,max(tmp)), na.color="light grey")(tmp), 
    border="light grey", lwd=0.5)
for (j in 1:177)
{
    nbrs = which(flux[,j,t[3]] > 0)
    if (!length(nbrs)) next
    wt = flux[nbrs,j,t[3]]
    nbrs = nbrs[which.max(wt)]
    wt = max(wt)
    m = colSums(sweep(site_coords[nbrs,,drop=FALSE], 1, wt, "*")) / sum(wt)
    arrows(
        0.85*m[1] + 0.15*site_coords[j,1], 
        0.85*m[2] + 0.15*site_coords[j,2], 
        0.85*site_coords[j,1] + 0.15*m[1], 
        0.85*site_coords[j,2] + 0.15*m[2], length=0.025,
        col=scales::alpha(1, 1000*sum(wt)), lwd=0.5)#lwd=sum(wt)*100)
}
title(sprintf("%s kya", times[t[3]]*25/1000), adj=0.1, cex.main=10/12)

dev.off()
