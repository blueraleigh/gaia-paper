#! /usr/local/bin/Rscript --vanilla

library(gaia)

ts = treeseq_load("../data/trees/chr18p.trees")
data = data.matrix(read.csv("../data/trees/sample-states.csv"))
cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))

mpr = treeseq_discrete_mpr(ts, data, cost.mat)

filename = "./results/mpr-chr18p.rds"

saveRDS(mpr, file=filename)
