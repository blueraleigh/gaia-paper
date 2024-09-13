library(sf)

landgrid = st_read("../../data/landgrid.gpkg")
land = st_transform(st_read("../../data/land.gpkg"), crs=st_crs(landgrid))

sites = st_centroid(landgrid$geom)

wrld = st_read("../../data/TM_WORLD_BORDERS-0.1.gpkg")
wrld = st_make_valid(wrld)

sites = st_transform(sites, st_crs(wrld))

SUBREGION = wrld$SUBREGION

site_map = wrld$SUBREGION[st_nearest_feature(sites, wrld)]
# map Taiwan to Eastern Asia
site_map[which(site_map == 0)] = 30L
# separate out Siberia east of Ural Mountains into a Northern Asia
site_map[site_map == 151L & (sapply(sites, "[", 1) > 59 | 
    sapply(sites, "[", 1) < -170)] = max(site_map) + 1L

sites = st_centroid(landgrid$geom)

site_coords = st_coordinates(sites)

# UN M49 regions
UNM49 = c(
`18`="#fe9929", # 5 southern africa
`17`="#ec7014", # 4 middle africa
`14`="#cc4c02", # 2 eastern africa
`11`="#993404", # 1 western africa
`15`="#662506", # 3 north africa

`39`="#4eb3d3",  # 9 southern europe
`151`="#2b8cbe", # 12 eastern europe
`154`="#0868ac", # 13 northern europe
`155`="#084081",  # 14 western europe

`30`="#99d8c9",  # 6 eastern asia
`34`="#66c2a4",  # 7 southern asia
`35`="#41ae76",  # 8 southeastern asia
`143`="#238b45", # 10 central asia
`145`="#006d2c", # 11 western asia
`156`="#00441b"  # 15 northern asia
)


ME = 145
EU = c(39,151,154,155)
AS = c(30,34,35,143,156)
AF = c(11,14,15,17,18)

ME_I = which(site_map %in% ME)
EU_I = which(site_map %in% EU)
AS_I = which(site_map %in% AS)
AF_I = which(site_map %in% AF)

ancestry = readRDS("../../analysis/results/ancestry-chr18p.rds")

Z = apply(ancestry, 2, function(a) {
    c(
    sum(a[AF_I]),
    sum(a[EU_I]),
    sum(a[AS_I]),
    sum(a[ME_I]))
})

times = attr(ancestry, "times")

plot(times, Z[1,], type='s', ylim=c(0,1), col="#a6611a", lwd=2, las=1,
    bty='l', ylab="Sample ancestry proportion", xlab="Generations ago")
lines(times, Z[2,], type='s', col="#dfc27d", lwd=2)
lines(times, Z[3,], type='s', col="#80cdc1", lwd=2)
lines(times, Z[4,], type='s', col="#018571", lwd=2)

legend("right", legend=c("Africa", "Europe", "Asia", "Middle East"),
    lty=1, col=c("#a6611a","#dfc27d","#80cdc1","#018571"), lwd=2, 
    bty="n")


dev.print(pdf, file="overall-ancestry.pdf")
