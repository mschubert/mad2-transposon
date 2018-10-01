library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')
gset = import('data/genesets')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('m', 'meta', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('i', 'cis', 'cis site RData', '../cis_analysis/poisson.RData'),
    opt('o', 'outfile', 'results RData', 'de.RData'),
    opt('p', 'plotfile', 'pdf', 'de.pdf'))

cis = io$load(args$cis)
cis_genes = cis$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name)
gene_copies = io$load(args$copies)
exprset = io$load(args$expr)
idx = io$load(args$meta) %>%
    mutate(type = ifelse(is.na(type), "unknown", type))
colnames(idx) = make.names(colnames(idx)) # fix "T-cell" in meta?
counts = exprset$counts
narray::intersect(gene_copies, counts, along=1)
narray::intersect(gene_copies, counts, idx$sample, along=2)

# vst w/ copy num corr
eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue+type) %>%
    DESeq2::estimateSizeFactors(normMatrix=gene_copies)
vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))

# fit tissue of origin and pan-aneuploidy
res = util$do_wald(eset, ~ tissue + type + aneuploidy, ex="tissue|aneuploidy")

# fit cancer type vs mean
cancer_type = function(term) {
    design(eset) = formula(paste("~ tissue + ", term))
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000) %>%
        extract_coef(term)
}
ats = c("T.cell", "Myeloid", "Other") # fix "T-cell" in meta?
res = c(res, sapply(ats, cancer_type, simplify=FALSE))

# fit aneuploidy within cancer types
aneup_tissue = function(type) {
    eset = eset[,eset[[type]] == 1]
    if (length(unique(eset$tissue)) == 1)
        design(eset) = formula(paste("~ aneuploidy"))
    else
        design(eset) = formula(paste("~ aneuploidy + tissue"))
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000) %>%
        extract_coef("aneuploidy")
}
res = c(res, setNames(lapply(ats, aneup_tissue), paste0(ats, ":aneuploidy")))

names(res) = sub("T\\.cell", "T-cell", names(res)) # fix "T-cell" in meta?

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(plot_pcs(idx, pca, 1, 2))
print(plot_pcs(idx, pca, 3, 4))
print(plot_pcs(idx, pca, 5, 6))
for (name in names(res)) {
    message(name)
    print(plot_volcano(res[[name]], cis_genes) + ggtitle(name))
}
dev.off()

save(eset, pca, cis, res, file=args$outfile)
