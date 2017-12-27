library(hdf5tenx)
library(BiocParallel)

fname <- "/home/mtmorgan/.ExperimentHub/1039"
group <- "mm10"

system.time({
    register(MulticoreParam(parallel::detectCores(), progressbar = TRUE))
    res <- margins(fname, group)
})

identical(sum(res[[2]][[1]] != 0), 1306127L)
identical(sum(as.numeric(res[[1]][[1]])), 2624828308)
all(mapply(function(x, y) {
    identical(sum(as.numeric(x)), sum(as.numeric(y)))
}, res[[1]], res[[2]]))
