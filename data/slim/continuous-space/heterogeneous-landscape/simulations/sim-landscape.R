library(fields)

# Simulate a Gaussian random field with an exponential covariance function,  
# range parameter = 2.0 and the domain is  [0,25] X [0,25] evaluating the 
# field at a 100x100 grid.  
grid = list(x=seq(0,25,,100), y=seq(0,25,,100)) 
obj = circulantEmbeddingSetup(grid, Covariance="Exponential", aRange=3)
for (i in 1:20) {
    look = circulantEmbedding(obj)
    write.table(
        sqrt(40)*look+30
        , file=sprintf("./landscapes/landscape-%d.csv", i)
        , row.names=FALSE
        , col.names=FALSE
        , sep=",")
}