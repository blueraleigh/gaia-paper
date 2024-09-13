#! /usr/local/bin/Rscript --vanilla

library(gaia)

NUMREPS = as.integer(commandArgs(trailingOnly=TRUE)[1])

ts = treeseq_load("../data/trees/chr18p.trees")
data = read.csv("../data/trees/sample-states.csv")
cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

mpr = readRDS("./results/mpr-chr18p.rds")

sample_sets = rep(1L, treeseq_num_samples(ts))
state_sets = 1:nrow(cost.mat)

times = seq(0, 20000, 100)

flux = matrix(0, length(state_sets), length(state_sets))

max.age = matrix(0, length(state_sets), length(state_sets))

flux.binned = array(0, c(length(state_sets), length(state_sets),
    length(times)-1))

for (i in 1:NUMREPS)
{
    tmp = treeseq_discrete_mpr_ancestry_flux(ts, mpr, cost.mat, 
        neighbor.mat, c(0, 20000), state_sets, sample_sets)

    flux = flux + tmp[,,1,1] / NUMREPS
    max.age = max.age + attr(tmp, "max.age") / NUMREPS

    tmp = treeseq_discrete_mpr_ancestry_flux(ts, mpr, 
        cost.mat, neighbor.mat, times, state_sets, sample_sets)

    flux.binned = flux.binned + tmp[,,1,] / NUMREPS
}

attr(flux, "max.age") = max.age

filename = "./results/flux-chr18p.rds"

saveRDS(flux, file=filename)

filename = "./results/flux-thru-time-chr18p.rds"

saveRDS(flux.binned, file=filename)
