library(hdf5tenx)

fname <- "/home/mtmorgan/.ExperimentHub/1039" # rle
group <- "mm10"

system.time({
    register(MulticoreParam(parallel::detectCores()))
    res <- margins_rle(fname, group)
})

str(res)
stopifnot(
    identical(sum(res$column$n != 0), 1306127L),
    identical(sum(as.numeric(res$column$n)), 2624828308),
    identical(sum(as.numeric(res$column$sum)), 6388703090),
    identical(sum(as.numeric(res$column$sumsq)), 270395442858),
    all(mapply(function(x, y) {
        identical(sum(as.numeric(x)), sum(as.numeric(y)))
    }, res$row, res$column))
)
