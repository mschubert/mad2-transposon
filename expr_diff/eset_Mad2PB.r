library(dplyr)
library(DESeq2)
io = import('io')
sys = import('sys')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('m', 'meta', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'eset_Mad2PB.RData'),
    opt('p', 'plotfile', 'pdf', 'eset_Mad2PB.pdf'))

gene_copies = io$load(args$copies)
exprset = io$load(args$expr)
idx = io$load(args$meta) %>%
    mutate(type = ifelse(is.na(type), "unknown", type))
colnames(idx) = sub("T-cell", "Tcell", colnames(idx))
counts = exprset$counts
narray::intersect(gene_copies, counts, along=1)
narray::intersect(gene_copies, counts, idx$sample, along=2)

# vst w/ copy num corr
eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue+type) %>%
    DESeq2::estimateSizeFactors(normMatrix=gene_copies)
vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

save(eset, file=args$outfile)
