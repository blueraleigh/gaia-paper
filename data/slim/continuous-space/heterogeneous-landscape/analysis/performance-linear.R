library(gaia)

locations = function(ts)
{
    ind = treeseq_individuals(ts)
    nodes = treeseq_nodes(ts)
    locs = t(apply(ind, 1, function(i) {
        c(i$individual_id, i$location[1:2])
    }))
    locs = cbind(
          nodes$node_id
        , nodes$is_sample
        , locs[match(nodes$individual_id, locs[, 1]), -1]
    )
    colnames(locs) = c("node_id","is_sample","x","y")
    locs
}

args = commandArgs(TRUE)

REP = args[1]
SUFFIX = args[2]

TREESEQ = file.path(getwd(), sprintf("../simulations/trees%s", SUFFIX), sprintf("tree-%s.trees", REP))

ts = treeseq_load(TREESEQ)
nodes = treeseq_nodes(ts)

idx = match(sample(unique(nodes$individual_id[nodes$is_sample == 1L]), 100L),
    nodes$individual_id)

ts2 = treeseq_simplify(ts, nodes$node_id[c(rbind(idx, idx+1L))])
nodes = treeseq_nodes(ts2)
edges = treeseq_edges(ts2)
inds = treeseq_individuals(ts2)
internal_nodes = which(nodes$is_sample != 1L)

locs = locations(ts2)
sample_locations = locs[locs[,2] == 1, c(1,3,4)]
ancestor_locations = locs[locs[,2] != 1,c(1,3,4)]
sample_centroid = colMeans(sample_locations[,2:3])

start = Sys.time()
mpr = treeseq_linear_mpr(ts2, sample_locations, TRUE, tol=1e-12)
stop = Sys.time()
elapsed = unclass(stop - start)[1]

map_x = treeseq_linear_mpr_minimize(mpr)

e = sqrt(rowSums((map_x[-(1:200),] - ancestor_locations[,2:3])^2)) / max(dist(sample_locations[,2:3]))

dist_from_sample_centroid0 = sqrt(
    rowSums(sweep(ancestor_locations[, 2:3], 2, sample_centroid)^2))

dist_from_sample_centroid = sqrt(
    rowSums(sweep(map_x[-(1:200), ], 2, sample_centroid)^2))

ans = data.frame(
    as.integer(REP),
    internal_nodes-1L,
    nodes[internal_nodes, "time"],
    e,
    dist_from_sample_centroid0,
    dist_from_sample_centroid

)

write.table(
    ans,
    file=sprintf("ancestor-estimates-linear%s.csv", SUFFIX),
    sep=",",
    col.names=FALSE,
    row.names=FALSE,
    append=TRUE
)

