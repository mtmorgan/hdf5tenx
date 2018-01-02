library(TENxBrainData)
library(ExperimentHub)
library(hdf5tenx)
library(data.table)

rle <- ExperimentHub()[["EH1039"]]
dense <- ExperimentHub()[["EH1040"]]

ncol <- 10000                           # 5000; 10000
                                        # n.b., 130 x 10k for full data
tenx <- TENxBrainData()[, seq_len(ncol)]

system.time({
    res1 <- colSums(as.matrix(counts(tenx)))
})                                      # 3.4s; 7.55s, 1Gb

options(DelayedArray.block.size=2e9)
system.time({
    res1a <- colSums(as.matrix(counts(tenx)))
})                                      # ; 6.9s, 1Gb

system.time({
    m <- h5read(dense, "counts", start = c(1, 1), count = dim(tenx))
    res2 <- colSums(m)
})                                      # 1.7s; 3.74s, 1Gb

system.time({
    res6 <- hdf5tenx:::margins_dense_slab(
        dense, "counts", ncol(tenx), 0, ncol(tenx)
    )
})                                      # ; 3.80, 743 Kb

system.time({
    start <- h5read(rle, "/mm10/indptr", start=1, count=ncol(tenx) + 1)
    idx <- cbind(
        h5read(rle, "/mm10/indices", start = 1, count = tail(start, 1)) + 1,
        rep(seq_len(length(start) - 1), diff(start))
    )
    m <- matrix(0L, nrow(tenx), ncol(tenx))
    m[idx] <- h5read(rle, "/mm10/data", start = 1, count = tail(start, 1))
    res3 <- colSums(m)
})                                      # 2.0; 4.3s, 1.3 Gb
    

system.time({
    count <- h5read(rle, "/mm10/indptr", start=ncol(tenx) + 1, count=1)
    dt <- data.table(
        row = h5read(rle, "/mm10/indices", start = 1, count = count) + 1,
        column = rep(seq_len(length(start) - 1), diff(start)),
        count = h5read(rle, "/mm10/data", start = 1, count = count)
    )
    res4 <- dt[, 
        list(n = .N, sum = sum(count), sumsq = sum(count * count)),
        keyby=column
    ]
})                                      # 0.76s; 1.48s, 310 Mb
                                        # 10k, sum only: 1.32s
                                        # 10k, row & col: 4.1s

system.time({
    offset <- h5read(rle, "/mm10/indptr", count = ncol(tenx) + 1)
    res5 <- hdf5tenx:::margins_rle_slab(rle, "mm10", 1e7, 0, ncol(tenx))
})                                      # 0.65; 1.29s, 743 Kb
print(object.size(res5) + object.size(res5$column[["0"]]), units="auto") 

