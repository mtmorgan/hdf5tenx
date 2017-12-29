library(hdf5tenx)

fname <- "/home/mtmorgan/.ExperimentHub/1040" # rectangular
group <- "counts"

system.time({
    register(MulticoreParam(parallel::detectCores(), progressbar = TRUE))
    res <- margins_dense(fname, group)
})

str(res)
stopifnot(
    identical(sum(res$column$n != 0), 1306127L),
    identical(sum(as.numeric(res$column$n)), 2624828308),
    all(mapply(function(x, y) {
        identical(sum(as.numeric(x)), sum(as.numeric(y)))
    }, res$row, res$column))
)
