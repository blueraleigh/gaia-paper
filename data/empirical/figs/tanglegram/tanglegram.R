library(sf)
library(gaia)

source("./bezier.R")

ts = treeseq_load("../../data/trees/chr18p.trees")
data = read.csv("../../data/trees/sample-states.csv")
mpr = readRDS("../../analysis/results/mpr-chr18p.rds")
h = read.csv("../../analysis/results/georef-arg.csv")
land = st_read("../../data/land.gpkg")
landgrid = st_read("../../data/landgrid.gpkg")

sites = st_centroid(landgrid$geom)

site_coords = st_coordinates(sites)
site_coords[6,] = st_coordinates(st_point_on_surface(landgrid$geom[6]))
site_coords[175,] = st_coordinates(st_point_on_surface(landgrid$geom[175]))

nodes = treeseq_nodes(ts)
edges = treeseq_edges(ts)

brks = classInt::classIntervals(nodes$time, style="fisher")
brks = classInt::classIntervals(h$time, style="fisher")
nb = length(brks$brks)
nb1 = nb - 1

pal = LSD::colorpalette("reds", max(range(classInt::findCols(brks))), rev=TRUE)

H = split(h, h$edge_id)

l = matrix(c(
    1,1,1,
    1,1,1,
    1,1,1,
    2,3,4,
    5,6,7), 5, 3, byrow=TRUE)


#jpeg("tanglegram.jpg", width=7.25,height=7.25,units="in",res=600)
tiff("tanglegram.tif", width=7.25,height=7.25,units="in",res=600)
layout(l)
par(mar=c(0,0,0,0))
plot.new()
plot.window(
    xlim=range(st_coordinates(land)[,1]),
    ylim=range(st_coordinates(land)[,2]), asp=1)
# work around macos sonoma plotting bug: https://stat.ethz.ch/pipermail/r-sig-mac/2024-February/014949.html
invisible(lapply(land$geom[[1]], function(L) polypath(sf:::p_bind(L), 
    lwd=0.5, rule="evenodd", col="light grey", border="light grey")))
invisible(sapply(1:length(H), function(i) {
    p = H[[i]]
    n = nrow(p)
    if (n > 2)
    {
        b = bezier(cbind(p$X[-1], p$Y[-1]), p$time[-1])(seq(0,1,.05))
        n = length(b$x)
        p = cbind(b$x[-n],b$x[-1],b$y[-n],b$y[-1])    
        segments(p[,1], p[,3],p[,2],p[,4],lwd=0.5,
            col=scales::alpha(pal[findInterval(b$time, brks$brks, all.inside=TRUE)],0.2))
    }
}))

plot(landgrid$geom, add=TRUE, lwd=0.15, border="grey")
plot(land$geom, add=TRUE, lwd=0.15)
points(site_coords, cex=0.25*sqrt(tabulate(data[,2], 177)), lwd=0.5,
    pch=21, bg=scales::alpha("white", 0.7))

# the rest of this really just makes the legend
g = st_make_grid(landgrid$geom, n=20, square=FALSE, flat_topped=FALSE)


bb = bezier(rbind(
    st_coordinates(st_centroid(g[300+8])),
    c(1382256.2, -2110351+1e6),
    st_coordinates(st_centroid(g[300-3*18-10]))
), c(0,0.5,1))(seq(0,1,,nb))
nn = length(bb$x)
pp = cbind(bb$x[-nn],bb$x[-1],bb$y[-nn]+2e5,bb$y[-1]+2e5)


segments(pp[,1], pp[,3],pp[,2],pp[,4],lwd=2)
arrows(pp[nb1,1], pp[nb1,3],pp[nb1,2],pp[nb1,4],lwd=2,
    length=0.05)
segments(pp[,1], pp[,3],pp[,2],pp[,4],lwd=1.8,col=rev(pal))
arrows(pp[nb1,1], pp[nb1,3],pp[nb1,2],pp[nb1,4],lwd=1.8,col=rev(pal)[nb1],
    length=0.05)

rect(
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,nb)[-1],
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-2e5,
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,nb)[-nb],
    st_coordinates(st_centroid(g[300+8]))[2],
    col=pal,
    border=pal
)
rect(
    st_coordinates(st_centroid(g[300-3*18-10]))[1],
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-2e5,
    st_coordinates(st_centroid(g[300+8]))[1],
    st_coordinates(st_centroid(g[300+8]))[2],lwd=0.5
)
segments(
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,5),
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-2e5,
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,5),
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-3e5,
    lwd=0.5
)
text(
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,5),
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-3e5,
    c("0","500","1000","1500","2000"), cex=0.6,
    offset=0.25, pos=1
)

text(
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,5)[3],
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-9e5,
    "Age of ancestral migration (kya)", cex=0.58
)


segments(
    st_coordinates(st_centroid(g[300-4*18+9]))[1],
    st_coordinates(st_centroid(g[300-4*18+9]))[2],
    st_coordinates(st_centroid(g[300-4*18+9]))[1] + 3.5e6,
    st_coordinates(st_centroid(g[300-4*18+9]))[2]
)
points(st_coordinates(st_centroid(g[300-4*18+9])), pch=21, bg="white")
points(
    seq(
        st_coordinates(st_centroid(g[300-4*18+9]))[1] + 3.5e6,
        st_coordinates(st_centroid(g[300-4*18+9]))[1] + 1e6, , 3),
    rep(st_coordinates(st_centroid(g[300-4*18+9]))[2], 3),
    pch=19
)

text(st_coordinates(st_centroid(g[300-4*18+9]))[1],
    st_coordinates(st_centroid(g[300-4*18+9]))[2]+3e5,
    "Path of average ancestor",
    cex=0.6, adj=0)
text(st_coordinates(st_centroid(g[300-4*18+9]))[1],
    st_coordinates(st_centroid(g[300-4*18+9]))[2], "0", pos=1, cex=0.6, adj=0)
text(
    seq(
        st_coordinates(st_centroid(g[300-4*18+9]))[1] + 3.5e6,
        st_coordinates(st_centroid(g[300-4*18+9]))[1] + 1e6, , 3),
    rep(st_coordinates(st_centroid(g[300-4*18+9]))[2], 3),
    c("2000 kya", "200", "10"),
    cex=0.6, adj=0, pos=1)

# here we make the subplots
times = brks$brks

show = c(1, 5, 289, 905, 2049, 2081)
for (sample_id in show)
{
    tss = treeseq_simplify(ts, sample_id-1L, node.map=TRUE,
        keep.unary=TRUE, keep.input.roots=TRUE, no.filter.nodes=TRUE)

    edges.subset = treeseq_edges(tss)
    edges.subset$edge_id = .Call(C_treeseq_fixup_edge_ids, tss@treeseq, 
        ts@treeseq, attr(tss, "node.map"))

    H2 = h[h$edge_id %in% edges.subset$edge_id,]
    H2$time.bin = cut(H2$time, times)

    H3 = split(H2, H2$edge_id)
    plot.new()
    plot.window(
        xlim=range(st_coordinates(land)[,1]),
        ylim=range(st_coordinates(land)[,2]), asp=1)
    invisible(lapply(land$geom[[1]], function(L) polypath(sf:::p_bind(L), 
        lwd=0.5, rule="evenodd", col="light grey", border="light grey")))
    
    invisible(sapply(1:length(H3), function(i) {
        p = H3[[i]]
        n = nrow(p)
        if (n > 2)
        {
            b = bezier(cbind(p$X[-1], p$Y[-1]), p$time[-1])(seq(0,1,.05))
            n = length(b$x)
            p = cbind(b$x[-n],b$x[-1],b$y[-n],b$y[-1])    
            segments(p[,1], p[,3],p[,2],p[,4],lwd=0.5,
                col=scales::alpha(pal[findInterval(b$time, brks$brks, all.inside=TRUE)],0.2))
        }
    }))
    plot(landgrid$geom, add=TRUE, lwd=0.15, border="grey")
    plot(land$geom, add=TRUE, lwd=0.15)
    centroid = t(sapply(split(H2, H2$time.bin), 
        function(df) colMeans(site_coords[df$state_id,,drop=FALSE])))
    idx = !apply(centroid, 1, anyNA)
    centroid = cbind(centroid[idx,], rowMeans(cbind(times[-nb],times[-1]))[idx])
    centroid = rbind(c(colMeans(site_coords[data[sample_id,2],,drop=FALSE]), 0), centroid)
    b = bezier(centroid[,1:2], centroid[,3])(seq(0,1,.01))
    lines(b$x, b$y)
    idx = sapply(c(400, 8000, 80000), 
        function(t) {
            i = findInterval(t, times)
            which(findInterval(b$time, times) == i)[1]
    })
    points(b$x[idx], b$y[idx], pch=19, cex=0.8)
    points(site_coords[data[sample_id,2],,drop=FALSE], pch=21, lwd=0.5,
        bg=scales::alpha("white", 1))
}

dev.off()
