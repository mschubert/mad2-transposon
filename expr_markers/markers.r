library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
    opt('m', 'meta', 'sample .RData', '../ploidy_compare/analysis_set.RData'),
    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

genes = c("Mad2l1", "Trp53", "Ets1", "Erg", "Myc", "Sox4", "Il7", "Il7r", "Kit", "Irf4", "Stat1", "Stat3",
          "Ebf1", "Cd19", "Ighm", "Ikzf1", "Tcf15", "Gata2", "Cd55b", "Pax5", "Tmem184a", "Dntt", "Igll1",
          "Vpreb1", "Rag1", "Rag2", "Xrcc6", "Cd34", "Cd3g", "Cd3e", "Klrb1c", "Tacr1", "Itgam", "Ly6a")

meta = io$load(args$meta)
counts = io$load(args$expr)$counts
rownames(counts) = idmap$gene(rownames(counts), to="mgi_symbol", dset="mmusculus_gene_ensembl")
counts = counts[rowSums(counts) > 0 & !is.na(rownames(counts)),]
narray::intersect(meta$sample, counts, along=2)

eset = DESeq2::DESeqDataSetFromMatrix(counts, meta, ~1) %>%
    DESeq2::estimateSizeFactors()
reads = DESeq2::counts(eset, normalized=TRUE)
expr = DESeq2::varianceStabilizingTransformation(eset)

exp1 = t(assay(expr[genes,]))
meta2 = meta %>% filter(type == "Other")
exp2 = t(assay(expr[genes, meta2$sample]))
pr = prcomp(exp1)
umap = uwot::umap(exp1)
pr2 = prcomp(exp2)
umap2 = uwot::umap(exp2)
colnames(umap) = colnames(umap2) = c("umap1", "umap2")
clust1 = igraph::cluster_louvain(scran::buildSNNGraph(t(exp1), k=5))
clust2 = igraph::cluster_louvain(scran::buildSNNGraph(t(exp2), k=5))
umap = cbind(meta, umap, clust=factor(clust1$membership))
umap2 = cbind(meta2, umap2, clust=factor(clust2$membership))

pdf("Rplots.pdf", 10, 8)
plt$pca(pr, aes(x=PC1, y=PC2), annot=meta, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(color=type, size=aneuploidy)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

ggplot(umap, aes(x=umap1, y=umap2)) +
    geom_point(aes(color=type, size=aneuploidy, shape=clust)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()

cbind(umap, t(reads[genes,])) %>%
    tidyr::gather("gene", "reads", -(aneup_src:clust)) %>%
    ggplot(aes(x=clust, y=log10(reads+1), group=clust)) +
        geom_point(aes(color=type, size=aneuploidy), alpha=0.3) +
        facet_wrap(~ gene, scales="free_y")

plt$pca(pr2, aes(x=PC1, y=PC2), annot=meta2, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(size=aneuploidy), color="#00ba38") +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

ggplot(umap2, aes(x=umap1, y=umap2)) +
    geom_point(aes(size=aneuploidy, shape=clust), color="#00ba38") +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()

cbind(umap2, t(reads[genes,meta2$sample])) %>%
    tidyr::gather("gene", "reads", -(aneup_src:clust)) %>%
    ggplot(aes(x=clust, y=log10(reads+1), group=clust)) +
        geom_point(aes(size=aneuploidy), color="#00ba38", alpha=0.3) +
        facet_wrap(~ gene, scales="free_y")
dev.off()
