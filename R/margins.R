#' @importFrom rhdf5 h5read
.iterator <- function(fname, group, bufsize = 1e7, n = Inf) {
    indname <- paste0(group, "/indptr")
    indptr <- h5read(fname, indname, bit64conversion = "double")
    indgrp <- floor(indptr / bufsize)
    offset <- match(unique(indgrp), indgrp)
    count <- as.integer(diff(indptr[c(offset, length(indptr))]))
    n <- min(n, length(count))
    i <- 0L
    function() {
        if (i == n)
            return(NULL)
        i <<- i + 1L
        list(fname = fname, group = group, offset=offset[i], count=count[i])
    }
}

#' @useDynLib hdf5tenx, .registration = TRUE
.fun <- function(iter, ...)
    margins_slab(iter$fname, iter$group, iter$offset - 1L, iter$count)

.reduce <- function(x, y)
    Map(Map, x, y, MoreArgs = list(f = `+`))

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

    bpiterate(.iterator(fname, group, bufsize, n), .fun, REDUCE=.reduce)
}
