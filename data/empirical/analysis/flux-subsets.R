#! /usr/local/bin/Rscript --vanilla

library(gaia)

NUMREPS = as.integer(commandArgs(trailingOnly=TRUE)[1])
num_subsets = 100L

cost.mat = data.matrix(read.csv("../data/landgrid-distmat.csv", row.names=1))
dimnames(cost.mat) = NULL
neighbor.mat = data.matrix(read.csv("../data/landgrid-adjmat.csv", row.names=1))
dimnames(neighbor.mat) = NULL

state_sets = 1:nrow(cost.mat)

times = seq(0, 20000, 100)

F = matrix(0, length(state_sets), length(state_sets))
Fv = matrix(0, length(state_sets), length(state_sets))

M = matrix(0, length(state_sets), length(state_sets))
Mv = matrix(0, length(state_sets), length(state_sets))

Fb = array(0, c(length(state_sets), length(state_sets),
    length(times)-1))
Fbv = array(0, c(length(state_sets), length(state_sets),
    length(times)-1))

for (j in 1:num_subsets)
{
    ts = treeseq_load(
        sprintf("../data/trees/subsets/chr18p-subset-%d.trees", j))

    mpr = readRDS(sprintf("./results/mpr-chr18p-subset-%d.rds", j))

    sample_sets = rep(1L, treeseq_num_samples(ts))

    flux = matrix(0, length(state_sets), length(state_sets))

    max.age = matrix(0, length(state_sets), length(state_sets))

    flux.binned = array(0, c(length(state_sets), length(state_sets),
        length(times)-1))

    for (i in 1:NUMREPS)
    {
        tmp = treeseq_discrete_mpr_ancestry_flux(ts, mpr, cost.mat, 
            neighbor.mat, c(0, 20000), state_sets, sample_sets)

        flux = flux + (tmp[,,1,1] - flux) / i
        max.age = max.age + (attr(tmp, "max.age") - max.age) / i

        tmp = treeseq_discrete_mpr_ancestry_flux(ts, mpr, 
            cost.mat, neighbor.mat, times, state_sets, sample_sets)

        flux.binned = flux.binned + (tmp[,,1,] - flux.binned) / i
    }

    Fm = F
    F = F + (flux - F) / j
    Fv = Fv + (flux - Fm)*(flux - F)

    Mm = M
    M = M + (max.age - M) / j
    Mv = Mv + (max.age - Mm)*(max.age - M)

    Fbm = Fb
    Fb = Fb + (flux.binned - Fb) / j
    Fbv = Fbv + (flux.binned - Fbm)*(flux.binned - Fb)

    rm(ts)
    rm(mpr)
    gc(full=TRUE)
}

Fv = Fv / (num_subsets - 1)
Mv = Mv / (num_subsets - 1)
Fbv = Fbv / (num_subsets - 1)

attr(F, "std.error") = sqrt(Fv) / sqrt(num_subsets)
attr(F, "max.age") = M
attr(F, "max.age.std.error") = sqrt(Mv) / sqrt(num_subsets)

attr(Fb, "std.error") = sqrt(Fbv) / sqrt(num_subsets)

filename = "./results/flux-subset-avg-chr18p.rds"

saveRDS(F, file=filename)

filename = "./results/flux-thru-time-subset-avg-chr18p.rds"

saveRDS(Fb, file=filename)
