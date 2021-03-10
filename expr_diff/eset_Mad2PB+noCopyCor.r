library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression rds', '../data/rnaseq/assemble.rds'),
    opt('m', 'meta', 'aneuploidy score', '../ploidy_compare/analysis_set.rds'),
    opt('o', 'outfile', 'results rds', 'eset_Mad2PB.rds'),
    opt('p', 'plotfile', 'pdf', 'eset_Mad2PB.pdf')
)

exprset = readRDS(args$expr)
idx = readRDS(args$meta) %>%
    mutate(type = ifelse(is.na(type), "unknown", type))
colnames(idx) = sub("T-cell", "Tcell", colnames(idx))
counts = exprset$counts
narray::intersect(counts, idx$sample, along=2)

# vst w/ copy num corr
eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue+type) %>%
    DESeq2::estimateSizeFactors()
rownames(eset) = idmap$gene(rownames(eset), from="ensembl_gene_id",
    to="external_gene_name", dset="mmusculus_gene_ensembl")
eset = eset[!is.na(rownames(eset)) & !duplicated(rownames(eset)),]
vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

saveRDS(list(eset=eset, vs=vs, pca=pca), file=args$outfile)
