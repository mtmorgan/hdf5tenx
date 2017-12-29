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
        list(begin = offset[i], end = offset[i + 1])
    }
}

.fun <- function(iter, ..., nrow, fname, group) {
    hdf5tenx:::margins_matrix(
        fname, group, nrow, iter$begin - 1L, iter$end - 1L
    )
}

.reduce <- function(x, y) {
    ## trying to avoid copying
    x$row$n <- x$row$n + y[[1]][[2]]
    x$row$sum <- x$row$sum + y[[1]][[3]]
    x$row$sumsq <- x$row$sumsq + y[[1]][[4]]

    idx <- y[[2]][[1]] + seq_along(y[[2]][[2]])
    x$column$n[idx] <- x$column$n[idx] + y[[2]][[2]]
    x$column$sum[idx] <- x$column$sum[idx] + y[[2]][[3]]
    x$column$sumsq[idx] <- x$column$sumsq[idx] + y[[2]][[4]]

    x
}

fname <- "/home/mtmorgan/.ExperimentHub/1040" # rectangular
group <- "counts"
dim <- hdf5tenx:::margins_dim(fname, group)

init <- list(
    row = list(
        n = integer(dim[1]), sum = numeric(dim[1]), sumsq = numeric(dim[1])
    ),
    column = list(
        n = integer(dim[2]), sum = numeric(dim[2]), sumsq = numeric(dim[2])
    )
)

register(MulticoreParam(parallel::detectCores(), progressbar = TRUE))
system.time({
    res <- bpiterate(
        .iterator(fname, group), .fun,
        nrow = dim[1], fname = fname, group = group,
        REDUCE = .reduce, init = init
    )
})

str(res)
stopifnot(
    identical(sum(res$column$n != 0), 1306127L),
    identical(sum(as.numeric(res$row$n)), 2624828308),
    all(mapply(function(x, y) {
        identical(sum(as.numeric(x)), sum(as.numeric(y)))
    }, res$row, res$column))
)
