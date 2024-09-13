library(sf)
library(gaia)

data = read.csv("./trees/sample-states.csv")
ts = treeseq_load("./trees/chr18p.trees")

landgrid = st_read("landgrid.gpkg")
land = st_read("land.gpkg")
sample.coords = st_as_sf(read.csv("./trees/sample-coords.csv"), 
    coords=c("X","Y"), crs=st_crs(land))

wrld = st_read("TM_WORLD_BORDERS-0.1.gpkg")
wrld = st_make_valid(wrld)

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

EUg = st_union(st_transform(
    wrld$geom[match(setdiff(wrld$NAME[wrld$SUBREGION %in% EU], "Russia"),
        wrld$NAME)], crs=st_crs(landgrid)))
AFg = st_union(st_transform(
    wrld$geom[wrld$SUBREGION %in% AF], 
    crs=st_crs(landgrid)))
MEg = st_union(st_transform(wrld$geom[wrld$SUBREGION %in% ME],
    crs=st_crs(landgrid)))
ASg = st_union(st_transform(
    wrld$geom[match(c(wrld$NAME[wrld$SUBREGION %in% AS], "Russia"),
        wrld$NAME)], crs=st_crs(landgrid)))

sample_map = integer(nrow(data))
sample_map[st_contains(MEg, sample.coords$geometry)[[1]]] = 4L
sample_map[st_contains(EUg, sample.coords$geometry)[[1]]] = 2L
sample_map[st_contains(ASg, sample.coords$geometry)[[1]]] = 3L
sample_map[st_contains(AFg, sample.coords$geometry)[[1]]] = 1L
sample_map[which(sample_map==0)] = sapply(which(sample_map==0), function(s) {
    st_nearest_feature(sample.coords$geometry[s], c(AFg,EUg,ASg,MEg))
})


nodes = treeseq_nodes(ts)
inds = treeseq_individuals(ts)
samples = t(sapply(inds[,4], function(b) jsonlite::fromJSON(rawToChar(b))))

sample_id = paste(samples[,"sample"], samples[,"sample_accession"], sep="_")

sample_region = data.frame(sample_id=sample_id[nodes$individual_id+1L], 
    sample_region=c("AF","EU","AS","ME")[sample_map])

write.csv(sample_region, file="./trees/sample-region.csv", row.names=FALSE)
