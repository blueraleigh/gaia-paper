#! /usr/local/bin/Rscript --vanilla

library(gaia)

NUMREPS = as.integer(commandArgs(trailingOnly=TRUE)[1])
num_subsets = 100L

cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

sample_region = read.csv("../data/trees/sample-region.csv")
sample_region[,2] = factor(sample_region[,2], levels=c("AF","EU","AS","ME"))

num_sample_sets = 4L
num_state_sets = nrow(cost.mat)

state_sets = 1:num_state_sets

times = seq(0, 20000, 100)

Q = array(0, c(num_state_sets, num_sample_sets, length(times)))
Qv = array(0, c(num_state_sets, num_sample_sets, length(times)))

for (i in 1:num_subsets)
{

    ts = treeseq_load(
        sprintf("../data/trees/subsets/chr18p-subset-%d.trees", i))
    data = read.csv(
        sprintf("../data/trees/subsets/sample-states-subset-%d.csv", i))

    nodes = treeseq_nodes(ts)
    inds = treeseq_individuals(ts)

    m = sapply(inds[,4], function(b) jsonlite::fromJSON(rawToChar(b)))
    sample_id = paste(m["sample", ], m["sample_accession",], sep="_")

    sample_map = sample_region[match(sample_id, sample_region[,1]), 2]
    
    sample_map = sample_map[nodes$individual_id[nodes$is_sample==1]+1]
    
    mpr = readRDS(sprintf("./results/mpr-chr18p-subset-%d.rds", i))

    ancestry = array(0, c(num_state_sets, num_sample_sets, length(times)))
    n = matrix(0, num_sample_sets, length(times))

    for (j in 1:NUMREPS)
    {
        tmp = treeseq_discrete_mpr_ancestry(ts, mpr, 
            cost.mat, neighbor.mat, times, state_sets, as.integer(sample_map))

        ancestry = ancestry + (tmp - ancestry) / j
    }

    Qm = Q
    Q = Q + (ancestry - Q) / i
    Qv = Qv + (ancestry - Qm)*(ancestry - Q)

    rm(ts)
    rm(mpr)
    gc(full=TRUE)
}

Qv = Qv / (num_subsets - 1)
attr(Qv, 'scale') = NULL

attr(Q, "std.error") = sqrt(Qv) / sqrt(num_subsets)
attr(Q, 'scale') = NULL

filename = "./results/ancestry-thru-time-subset-avg-chr18p.rds"

saveRDS(Q, file=filename)

