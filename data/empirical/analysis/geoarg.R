library(gaia)
library(sf)

ts = treeseq_load("../data/trees/chr18p.trees")
data = read.csv("../data/trees/sample-states.csv")
cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

mpr = readRDS("./results/mpr-chr18p.rds")

num_states = nrow(cost.mat)

landgrid = st_read("../data/landgrid.gpkg")

sites = st_coordinates(st_centroid(landgrid))
sites[6,] = st_coordinates(st_point_on_surface(landgrid$geom[6]))
sites[175,] = st_coordinates(st_point_on_surface(landgrid$geom[175]))

edges = treeseq_edges(ts)

# sample a migration history for each edge
h = treeseq_discrete_mpr_edge_history(ts, mpr, cost.mat, neighbor.mat)
node.state = attr(h, 'node.state')

# all nodes get mapped to hex cell centers for pretty plotting later 
XY.node = sites[node.state,]

# assign each lineage a set of random coordinates within hex cells on its
# migration path. these will be used as control points for plotting
xy.edge = st_coordinates(st_sample(landgrid, tabulate(h[,2], num_states)))
stopifnot(nrow(xy.edge) == nrow(h))
XY.edge = matrix(0, nrow(xy.edge), 2)
XY.edge[order(h[,2]), ] = xy.edge

# map the endpoint locations of each edge back to hex cell centers
path_start_index = head(attr(h, "path.offset"), -1)
path_end_index = head(attr(h, "path.offset"), -1) + 
    diff(attr(h, "path.offset")) - 1L
XY.edge[path_start_index,] = XY.node[edges$child_id + 1L, ]
XY.edge[path_start_index+1,] = XY.node[edges$child_id + 1L, ]
XY.edge[path_end_index,] = XY.node[edges$parent_id + 1L, ]

h$X = XY.edge[,1]
h$Y = XY.edge[,2]


filename = "./results/georef-arg.csv"

write.csv(h, file=filename, row.names=FALSE)
