library(hdf5tenx)

.iterator <- function(fname, group, bufsize = 6000, n = Inf) {
    dim <- hdf5tenx:::margins_dim(fname, group)
    colgrp <- ceiling(seq_len(dim[[2]]) / bufsize)
    offset <- c(match(unique(colgrp), colgrp), length(colgrp) + 1L)
    n <- min(n, length(offset) - 1L)
    i <- 0L
    function() {
        if (i == n)
            return(NULL)
        i <<- i + 1L
        list(
            fname = fname, group = group,
            nrow = dim[1], begin = offset[i], end = offset[i + 1]
        )
    }
}

.fun <- function(iter, ...) {
    hdf5tenx:::margins_matrix(
        iter$fname, iter$group, iter$nrow, iter$begin - 1L, iter$end - 1L
    )
}

.reduce <- hdf5tenx:::.reduce

.final <- hdf5tenx:::.final

fname <- "/home/mtmorgan/.ExperimentHub/1040" # rectangular
group <- "counts"

register(MulticoreParam(parallel::detectCores(), progressbar = TRUE))
system.time({
    res <- bpiterate(
        .iterator(fname, group), .fun, REDUCE = .reduce
    )

    .final(res)
})

str(res)
stopifnot(
    identical(sum(res$column$n != 0), 1306127L),
    identical(sum(as.numeric(res$row$n)), 2624828308),
    all(mapply(function(x, y) {
        identical(sum(as.numeric(x)), sum(as.numeric(y)))
    }, res$row, res$column))
)
