source("data.R")

set.seed(1123581321L)

num_subsets = 100L

D = data_full(18L)

treeseq_write(D$ts, "./trees/chr18p.trees")
write.csv(D$data, file="./trees/sample-states.csv", row.names=FALSE)
write.csv(D$sample.coords, file="./trees/sample-coords.csv", row.names=FALSE)

for (i in 1:num_subsets)
{
    D2 = data_subsample(D)
    treeseq_write(D2$ts, sprintf("./trees/subsets/chr18p-subset-%d.trees", i))
    write.csv(D2$data, 
        file=sprintf("./trees/subsets/sample-states-subset-%d.csv", i), 
        row.names=FALSE)
    write.csv(D2$sample.coords, 
        file=sprintf("./trees/subsets/sample-coords-subset-%d.csv", i), 
        row.names=FALSE)
}
