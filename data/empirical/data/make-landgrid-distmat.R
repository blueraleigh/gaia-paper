library(sf)

setwd("~/Dropbox/research/postdoc/bradburd/manuscripts/gaia/v2-science-revision/hgdp/data")

source("proj.R")

x = st_read("landgrid.gpkg")

#D = st_distance(st_centroid(x))
D = st_distance(st_centroid(st_transform(x, GRS80)))

z = data.matrix(read.csv("landgrid-adjmat.csv", row.names=1))
dimnames(z) = NULL

G = igraph::graph_from_adjacency_matrix(data.matrix(z), mode="undirected", 
    weighted=TRUE)


num_edges = length(igraph::E(G))
el = igraph::as_edgelist(G)


for (i in 1:num_edges) {
    e = el[i,]
    from = e[1]
    to = e[2]
    igraph::E(G)[i]$weight = D[from, to]
}


igraph::E(G)$weight = igraph::E(G)$weight / 1000000

d = igraph::distances(G)

#write.csv(d, file="landgrid-distmat-km-proj.csv")
#write.csv(d^2, file="landgrid-distmat-squared-km-proj.csv")
write.csv(d, file="landgrid-distmat-km-geodesic.csv")
write.csv(d^2, file="landgrid-distmat-squared-km-geodesic.csv")

a = d
a[which(z == 0, arr.ind=TRUE)] = 0
write.csv(a, file="landgrid-adjmat-km-geodesic.csv")

plot(x$geom)

#centers = st_coordinates(st_centroid(x))
#for (i in 1:177) {
#    for (j in i:177) {
#        if (z[i,j]==1)
#        {
#            segments(centers[i,1],centers[i,2],centers[j,1],centers[j,2],
#               col=scales::col_numeric(LSD::colorpalette("spectral", rev=TRUE), 
#               range(d[which(z==1,arr.ind=TRUE)]))(d[i,j]))
#        }
#    }
#}

