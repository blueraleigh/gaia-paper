library(sf)
library(gaia)

data = read.csv("../../data/trees/sample-states.csv")
#Z = readRDS("../../analysis/results/ancestry-thru-time-subset-avg-chr18p.rds")
Z = readRDS("../../analysis/results/ancestry-thru-time-chr18p.rds")

landgrid = st_read("../../data/landgrid.gpkg")
land = st_read("../../data/land.gpkg")
sample.coords = st_as_sf(read.csv("../../data/trees/sample-coords.csv"), 
    coords=c("X","Y"), crs=st_crs(land))

sample_region = read.csv("../../data/trees/sample-region.csv")

sample_map = as.integer(factor(sample_region[,2], levels=c("AF","EU","AS","ME")))

sites = st_centroid(landgrid$geom)

site_coords = st_coordinates(sites)


par(mfcol=c(5,4), mar=c(0,0,0,0))
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,1,1]))
plot(land$geom, col="light grey", border="light grey", lwd=0.5)
points(site_coords, cex=0.1*sqrt(tabulate(data[sample_map==1,2], 177)), lwd=0.5,
    pch=21, bg=scales::alpha("white", 0.9))
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,1,5]))
plot(landgrid$geom, col=pal(Z[,1,5]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,1,41]))
plot(landgrid$geom, col=pal(Z[,1,41]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,1,81]))
plot(landgrid$geom, col=pal(Z[,1,81]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,1,201]))
plot(landgrid$geom, col=pal(Z[,1,201]), border="light grey", lwd=0.5)

pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,2,1]))
plot(land$geom, col="light grey", border="light grey", lwd=0.5)
points(site_coords, cex=0.1*sqrt(tabulate(data[sample_map==2,2], 177)), lwd=0.5,
    pch=21, bg=scales::alpha("white", 0.9))
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,2,5]))
plot(landgrid$geom, col=pal(Z[,2,5]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,2,41]))
plot(landgrid$geom, col=pal(Z[,2,41]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,2,81]))
plot(landgrid$geom, col=pal(Z[,2,81]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,2,201]))
plot(landgrid$geom, col=pal(Z[,2,201]), border="light grey", lwd=0.5)


pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,3,1]))
plot(land$geom, col="light grey", border="light grey", lwd=0.5)
points(site_coords, cex=0.1*sqrt(tabulate(data[sample_map==3,2], 177)), lwd=0.5,
    pch=21, bg=scales::alpha("white", 0.9))
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,3,5]))
plot(landgrid$geom, col=pal(Z[,3,5]), border="light grey", lwd=0.5)
g = st_make_grid(landgrid$geom, square=FALSE, flat_topped=TRUE,
    n=20)
plot(g[rev(c(107,131,155,179,203,227))], add=TRUE, border="light grey",
    col=scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
        c(0,1), na.color="light grey")(c(0.05, 0.25, 0.45, 0.65, 0.85,1)),
        lwd=0.5)
text(
    st_coordinates(st_centroid(g[227]))[1],
    st_coordinates(st_centroid(g[227]))[2],
    "Low ancestry", pos=4, cex=0.8, xpd=NA)
text(
    st_coordinates(st_centroid(g[107]))[1],
    st_coordinates(st_centroid(g[107]))[2],
    "High ancestry", pos=4, cex=0.8, xpd=NA)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,3,41]))
plot(landgrid$geom, col=pal(Z[,3,41]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,3,81]))
plot(landgrid$geom, col=pal(Z[,3,81]), border="light grey", lwd=0.5)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,3,201]))
plot(landgrid$geom, col=pal(Z[,3,201]), border="light grey", lwd=0.5)


pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,4,1]))
plot(land$geom, col="light grey", border="light grey", lwd=0.5)
legend("bottomright", legend="contemporary samples", cex=0.8, bty="n", inset=0.1)
points(site_coords, cex=0.1*sqrt(tabulate(data[sample_map==4,2], 177)), lwd=0.5,
    pch=21, bg=scales::alpha("white", 0.9))
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,4,5]))
plot(landgrid$geom, col=pal(Z[,4,5]), border="light grey", lwd=0.5)
legend("bottomright", legend="10 kya", cex=0.8, bty="n", inset=0.1)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,4,41]))
plot(landgrid$geom, col=pal(Z[,4,41]), border="light grey", lwd=0.5)
legend("bottomright", legend="100 kya", cex=0.8, bty="n", inset=0.1)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,4,81]))
plot(landgrid$geom, col=pal(Z[,4,81]), border="light grey", lwd=0.5)
legend("bottomright", legend="200 kya", cex=0.8, bty="n", inset=0.1)
pal = scales::col_numeric(c("light grey", LSD::colorpalette("reds", rev=TRUE)),
    range(Z[,4,201]))
plot(landgrid$geom, col=pal(Z[,4,201]), border="light grey", lwd=0.5)
legend("bottomright", legend="500 kya", cex=0.8, bty="n", inset=0.1)

#dev.print(pdf, file="ancestry-thru-time-subset-avg.pdf")
dev.print(pdf, file="ancestry-thru-time.pdf")
