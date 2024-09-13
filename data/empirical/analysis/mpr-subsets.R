#! /usr/local/bin/Rscript --vanilla

library(gaia)

cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
num_subsets = 100L

for (i in 1:num_subsets)
{
    ts = treeseq_load(
        sprintf("../data/trees/subsets/chr18p-subset-%d.trees", i))
    data = data.matrix(read.csv(
        sprintf("../data/trees/subsets/sample-states-subset-%d.csv", i)))

    mpr = treeseq_discrete_mpr(ts, data, cost.mat)

    filename = sprintf("./results/mpr-chr18p-subset-%d.rds", i)
    
    saveRDS(mpr, file=filename)
    rm(ts)
    rm(mpr)
    gc(full=TRUE)
}
