#! /usr/local/bin/Rscript --vanilla

library(gaia)

source('true-flux.R')

args = commandArgs(trailingOnly=TRUE)
REP = args[1]
MAP = args[2]

data = read.csv("../../../hgdp/data/trees/sample-states.csv")
pops_to_sample = sort(unique(data[,2])-1L)

cost.mat = data.matrix(read.csv("../../../hgdp/data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../../../hgdp/data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

ts = treeseq_load(sprintf("../simulations/trees/tree-ooa-%s-%s.trees", MAP, REP))
nodes = treeseq_nodes(ts)
sample_nodes = nodes[nodes$is_sample==1, ]
sample_nodes_to_sample = sample_nodes[sample_nodes$population_id %in% pops_to_sample,]
indivs = unname(tapply(sample_nodes_to_sample$individual_id, sample_nodes_to_sample$population_id,
    sample, 1))
to_keep = sample_nodes_to_sample$node_id[sample_nodes_to_sample$individual_id %in% indivs]

ts2 = treeseq_simplify(ts, to_keep, filter.populations=FALSE)

nodes = treeseq_nodes(ts2)
sample_nodes = nodes[nodes$is_sample==1, ]

sample_locations = cbind(node_id=sample_nodes[,1], state_id=sample_nodes[,4]+1L)

mpr = treeseq_discrete_mpr(ts2, sample_locations, cost.mat)
estimated_node_states = treeseq_discrete_mpr_minimize(mpr)

write.csv(data.frame(node_time=nodes$time, node_state=nodes$population_id+1L, 
    estimated_node_state=estimated_node_states), 
    file=sprintf("./mpr2/ooa-%s-mpr-%s.csv", MAP, REP), row.names=FALSE)


sample_sets = rep(1L, treeseq_num_samples(ts2))
state_sets = 1:nrow(cost.mat)

flux = treeseq_discrete_mpr_ancestry_flux(ts2, mpr, cost.mat, 
    neighbor.mat, c(0, 20000), state_sets, sample_sets)
attr(flux, "sites") = unique(sample_locations[,2])

ts3 = treeseq_simplify(ts, to_keep, filter.populations=FALSE, keep.unary=TRUE)
nodes0 = treeseq_nodes(ts3)

flux0 = true_flux(nodes0$population_id, ts3, c(0,20000), cost.mat, neighbor.mat,
    sample_sets, state_sets)

f = data.frame(true_flux=c(flux0), estimated_flux=c(flux))

write.csv(f, file=sprintf("./flux2/ooa-%s-flux-%s.csv", MAP, REP))

ooa_flux = data.frame(
    scenario=MAP,
    rep=REP,
    flux_north=sum(flux[38,69:70,1,1]), 
    flux_south=flux[67,68,1,1])

write_col_names = (REP == "1") && (MAP == "north")
write.table(ooa_flux, file='./flux2/ooa-flux.csv', sep=",", row.names=FALSE,
    col.names=ifelse(write_col_names, TRUE, FALSE),
    append=TRUE)


