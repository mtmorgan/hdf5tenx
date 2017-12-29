.iterator <- function(fname, group, bufsize = 1e7, n = Inf) {
    indptr <- indptr(fname, group)
    indgrp <- floor(indptr / bufsize)
    offset <- c(match(unique(indgrp), indgrp), length(indptr))
    n <- min(n, length(offset) - 1L)
    i <- 0L
    function() {
        if (i == n)
            return(NULL)
        i <<- i + 1L
        list(
            fname = fname, group = group,
            indptr = indptr, begin = offset[i], end = offset[i + 1]
        )
    }
}

#' @useDynLib hdf5tenx, .registration = TRUE
#' @importFrom Rcpp evalCpp
.fun <- function(iter, ...)
    margins_slab(
        iter$fname, iter$group, iter$indptr, iter$begin - 1L, iter$end - 1L
    )

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

#' Calculate row and column counts, sums, and sums-of-squares
#'
#' @param fname character(1) HDF5 file name to 10xGenomics-format file.
#'
#' @param group character(1) group name (e.g., `mm10`) containing counts.
#'
#' @param bufsize integer(1) approximate maximum number of counts to
#'     read per iteration. The default is reasonably memory- and
#'     time-efficient.
#'
#' @param n numeric(1) number of iterations to perform; default: all.
#'
#' @return a `list()` with `row` and `column` elements. Each element
#'     is itself a list with `n`, `sum`, and `sumsq` (sum-of-squares)
#'     of counts in the corresponding row or column.
#'
#' @importFrom BiocParallel bpiterate
#'
#' @export
margins <- function(fname, group, bufsize = 1e7, n = Inf) {
    stopifnot(
        is.character(fname), length(fname) == 1L, !is.na(fname),
        file.exists(fname),
        is.character(group), length(group) == 1L, !is.na(group),
        is.numeric(bufsize), length(bufsize) == 1L, bufsize > 0,
        is.numeric(n), length(bufsize) == 1L, n > 0
    )
    dim <- margins_dim(fname, group)
    init <- list(
        row = list(
            n = integer(dim[1]), sum = numeric(dim[1]), sumsq = numeric(dim[1])
        ),
        column = list(
            n = integer(dim[2]), sum = numeric(dim[2]), sumsq = numeric(dim[2])
        )
    )
    bpiterate(
        .iterator(fname, group, bufsize, n), .fun, REDUCE=.reduce, init = init
    )
}
