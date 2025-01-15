library(sf)

source("./arrow.R")

data = read.csv("../../data/trees/sample-states.csv")
landgrid = st_read("../../data/landgrid.gpkg")
land = st_read("../../data/land.gpkg")
site_coords = st_coordinates(st_centroid(landgrid))
site_coords = st_coordinates(st_point_on_surface(landgrid))

flux = readRDS("../../analysis/results/flux-chr18p.rds")
M = attr(flux, "max.age")

colfn = scales::col_numeric(LSD::colorpalette("reds", rev=TRUE), 
    range(M[M > 0]))

pdf(file="overall-flux.pdf", width=7, height=7, colormodel='cmyk')

par(mar=c(0,0,0,0), ps=9)
plot(landgrid$geom, border="light grey", lwd=0.5)
plot(land$geom, add=TRUE, lwd=0.5)

s = unique(data[,2])

for (ii in s)
{
    prev = -1L
    state_id = ii
    repeat {
        i = which.max(flux[,state_id])
        from.x = site_coords[i,1]
        from.y = site_coords[i,2]
        to.x = site_coords[state_id,1]
        to.y = site_coords[state_id,2]
        arrow(
            0.9*from.x + 0.1*to.x,
            0.9*from.y + 0.1*to.y,
            0.9*to.x + 0.1*from.x,
            0.9*to.y + 0.1*from.y,
            col=colfn(M[i,state_id]),
            awd=1e5,
            length=1e5,
            f=0.01,
            border=1,
            lwd=0.25
        )
        prev = state_id
        state_id = i
        if (which.max(flux[,state_id]) == prev)
            break
    }
}

points(
    site_coords[sort(unique(data[,2])),],
    cex=0.25*sqrt(table(data[,2])),
    lwd=0.5,
    pch=21,
    bg=scales::alpha("white", 0.5)
)


g = st_make_grid(landgrid, n=20, square=FALSE, flat_topped=FALSE)


plot(g[300-3*18], add=TRUE, border="light grey")

arrow(
    st_coordinates(st_centroid(g[300-4*18+9]))[1],
    st_coordinates(st_centroid(g[300-4*18+9]))[2],
    st_coordinates(st_centroid(g[300-3*18]))[1],
    st_coordinates(st_centroid(g[300-3*18]))[2],
    awd=1e5, length=1e5, lwd=0.25, f=0.01, col="white"
)
text(
    st_coordinates(st_centroid(g[300-4*18+9]))[1],
    st_coordinates(st_centroid(g[300-4*18+9]))[2],
    "Direction of\ngreatest migration", cex=1, pos=2)
text(
    st_coordinates(st_centroid(g[300]))[1]+9e5,
    st_coordinates(st_centroid(g[300]))[2],
    "Age of earliest migration\nancestral to sample (kya)", cex=.8,
    adj=1
)
rect(
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,101)[-1],
    st_coordinates(st_centroid(g[300-3*18-10]))[2]-2e5,
    seq(
        st_coordinates(st_centroid(g[300-3*18-10]))[1],
        st_coordinates(st_centroid(g[300+8]))[1],,101)[-101],
    st_coordinates(st_centroid(g[300+8]))[2],
    col=colfn(seq(min(M[M>0]), max(M[M>0]),,100)),
    border=NA
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
    c("0","125","250","375","500"), cex=1,
    offset=0.25, pos=1
)


# https://stat.ethz.ch/pipermail/r-help/2003-July/035796.html
# for converting cex to plot ("usr") coordinates
scal1 = (par("cin")[2] / par("pin")[1]) * 
    (par("usr")[2] - par("usr")[1]) * 0.5 * 0.375
scal2 = (par("cin")[2] / par("pin")[2]) * 
    (par("usr")[4] - par("usr")[3]) * 0.5 * 0.375

pt.cex=0.25*sqrt(c(355, 205, 55))
pt.xrad = pt.cex*scal1
pt.yrad = pt.cex*scal2


px = st_coordinates(st_centroid(g[300-3*18]))[1]
py = rep(st_coordinates(st_centroid(g[300-3*18]))[2]-2.2e6, 3)
px = c(px, px+pt.xrad[1] - pt.xrad[2], px+pt.xrad[1] - pt.xrad[3])

points(px, py, cex=pt.cex, lwd=0.5)
text(
    px[1]+pt.xrad[1],
    py[1],
    "Sample size", pos=4, cex=1)

segments(
    px[1]+pt.xrad[1],
    py[1]-pt.yrad[1]-1e5,
    px[1]-pt.xrad[1],
    py[1]-pt.yrad[1]-1e5,
    lwd=0.5)
segments(
    px[1]+pt.xrad[1]-0.25*sqrt(c(0,55,205,355))*scal1*2,
    py[1]-pt.yrad[1]-1e5,
    px[1]+pt.xrad[1]-0.25*sqrt(c(0,55,205,355))*scal1*2,
    py[1]-pt.yrad[1]-2e5, lwd=0.5
)
text(
    px[1]+pt.xrad[1]-0.25*sqrt(c(0,55,205,355))*scal1*2,
    py[1]-pt.yrad[1]*c(1.08,1.25,1.4,1.4),
    c(0,50,200,350),
    pos=1, cex=.8, srt=90
)

dev.off()
