library(dplyr)
io = import('io')

tps = io$load("../../data/rnaseq/assemble.RData")$counts
p53 = io$load("../../data/rnaseq/Mad2+p53_batch1.RData")$counts

mat = narray::stack(tps, p53, along=2)

expr = DESeq2::DESeqDataSetFromMatrix(mat, colData=data.frame(name=colnames(mat)), ~1) %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::counts(normalized=TRUE)

saveRDS(expr, file="expr.rds")
