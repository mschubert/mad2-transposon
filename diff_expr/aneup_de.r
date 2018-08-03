library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'aneup_de.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_de.pdf'))

gene_copies = io$load(args$copies)
exprset = io$load(args$expr)
idx = exprset$idx %>%
    mutate(sample = paste0(hist_nr, tissue)) %>%
    left_join(io$load(args$aneup)) %>%
    mutate(type = ifelse(is.na(type), "NA", tissue), #TODO: add annotations
           tissue = factor(tissue), type=factor(type))
counts = exprset$counts
narray::intersect(gene_copies, counts, along=1)
narray::intersect(gene_copies, counts, idx$sample, along=2)

# vst w/ copy num corr
eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue) %>%
    DESeq2::estimateSizeFactors(normMatrix=gene_copies) %>%
    DESeq2::estimateDispersions()

# pca plts
vs = DESeq2::getVarianceStabilizedData(eset)
pca = prcomp(t(vs), scale=F) #TRUE)
ggplot(cbind(idx, pca$x), aes(x=PC1, y=PC2, color=type, shape=tissue)) +
    geom_point(aes(size=aneup)) +
    geom_text_repel(aes(label=sample), color="black") +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
         title = "PCA plot")

design(eset) = ~ tissue + type + aneup
res1 = DESeq2::DESeq(eset) #, test="LRT", full=~tissue+type, reduced=~tissue)
resultsNames(res1)
#tab = lfcShrink(res1, coef="tissue_t_vs_s") %>% #, type="apeglm") %>%
tab = lfcShrink(res1, coef="aneup") %>% #, type="apeglm") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    tbl_df() %>%
    arrange(pvalue)

# res2 = DESeq2::DESeq(eset, test="LRT", reduced=~tissue + type) # higher pvals overall

pdf(args$plotfile)

ggplot(cbind(idx, pca$x), aes(x=PC3, y=PC4, color=type, shape=tissue)) +
    geom_point(aes(size=aneup)) +
    geom_text_repel(aes(label=sample), color="black") +
    labs(x = sprintf("PC3 (%.1f%%)", summary(pca)$importance[2,3]*100),
         y = sprintf("PC4 (%.1f%%)", summary(pca)$importance[2,4]*100),
         title = "PCA plot")

ggplot(cbind(idx, pca$x), aes(x=PC5, y=PC6, color=type, shape=tissue)) +
    geom_point(aes(size=aneup)) +
    geom_text_repel(aes(label=sample), color="black") +
    labs(x = sprintf("PC5 (%.1f%%)", summary(pca)$importance[2,5]*100),
         y = sprintf("PC6 (%.1f%%)", summary(pca)$importance[2,6]*100),
         title = "PCA plot")
dev.off()

# select number of PCs to regress out

# do DEseq on counts w/ normmatrix
