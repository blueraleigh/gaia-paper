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

SIGMA = args[1]
REP = args[2]

TREESEQ = file.path(getwd(), "../simulations/trees", sprintf("tree-S%s-R%s.trees", SIGMA, REP))

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

mpr = treeseq_linear_mpr(ts2, sample_locations, TRUE)

mean_map_rate = sqrt(mpr[[1]] / 2)

effective_rate = sqrt(mean(apply(data.matrix(edges), 1L, function(e) {
    parent = nodes[e[4]+1, 5]
    child = nodes[e[5]+1, 5]
    len = nodes[e[4]+1, 3] - nodes[e[5]+1, 3]
    from = inds[parent+1, ]$location[1:2]
    to = inds[child+1, ]$location[1:2]
    sum(abs(from - to)) / len
})) / 2)

ans = data.frame(
    as.numeric(SIGMA),
    as.integer(REP),
    mean_map_rate,
    effective_rate
)
write.table(
    ans,
    file="rate-estimates.csv",
    sep=",",
    col.names=FALSE,
    row.names=FALSE,
    append=TRUE
)


map_x = treeseq_linear_mpr_minimize(mpr)

e = rowSums(abs(map_x[-(1:200),] - ancestor_locations[,2:3])) / max(dist(sample_locations[,2:3], method="manhattan"))

ans = data.frame(
    as.numeric(SIGMA),
    as.integer(REP),
    internal_nodes-1L,
    nodes[internal_nodes, "time"],
    e
)

write.table(
    ans,
    file="ancestor-estimates.csv",
    sep=",",
    col.names=FALSE,
    row.names=FALSE,
    append=TRUE
)
