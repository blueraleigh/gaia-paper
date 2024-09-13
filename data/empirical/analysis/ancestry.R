#! /usr/local/bin/Rscript --vanilla

library(gaia)

NUMREPS = as.integer(commandArgs(trailingOnly=TRUE)[1])

ts = treeseq_load("../data/trees/chr18p.trees")
data = read.csv("../data/trees/sample-states.csv")
cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

sample_region = read.csv("../data/trees/sample-region.csv")
sample_region[,2] = factor(sample_region[,2], levels=c("AF","EU","AS","ME"))

num_sample_sets = 4L
num_state_sets = nrow(cost.mat)

mpr = readRDS("./results/mpr-chr18p.rds")

num_sample_sets = 4L
num_state_sets = nrow(cost.mat)

sample_sets = as.integer(sample_region[,2])
state_sets = 1:num_state_sets

times = seq(0, 20000, 100)

ancestry = array(0, c(num_state_sets, num_sample_sets, length(times)))
n = matrix(0, num_sample_sets, length(times))

for (i in 1:NUMREPS)
{

    tmp = treeseq_discrete_mpr_ancestry(ts, mpr, 
        cost.mat, neighbor.mat, times, state_sets, sample_sets)

    ancestry = ancestry + (tmp - ancestry) / i
}

filename = "./results/ancestry-thru-time-chr18p.rds"

saveRDS(ancestry, file=filename)

sample_sets = rep(1L, treeseq_num_samples(ts))

times = seq(0, 82000, 500)

ancestry = matrix(0, length(state_sets), length(times))
n = numeric(length(times))

for (i in 1:NUMREPS)
{

    tmp = treeseq_discrete_mpr_ancestry(ts, mpr, 
        cost.mat, neighbor.mat, times, state_sets, sample_sets)

    n = n + (attr(tmp, "scale")[1,] - n) / i

    ancestry = ancestry + (tmp[,1,] - ancestry) / i
}

ancestry = ancestry[,-165]
attr(ancestry, "times") = times[-165]
attr(ancestry, "scale") = n[-165]

saveRDS(ancestry, file="./results/ancestry-chr18p.rds")
