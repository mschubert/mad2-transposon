library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression rds', '../data/rnaseq/assemble.rds'),
    opt('m', 'meta', 'sample rds', '../ploidy_compare/analysis_set.rds'),
    opt('o', 'outfile', 'rds', 'markers.rds'),
    opt('p', 'plotfile', 'pdf to save plot to', 'markers.pdf')
)

genes = c("Mad2l1", "Trp53", "Ets1", "Erg", "Myc", "Sox4", "Il7", "Il7r", "Kit", "Irf4", "Stat1", "Stat3", "Bcl6",
          "Ebf1", "Cd19", "Ighm", "Ikzf1", "Tcf15", "Gata2", "Cd55b", "Pax5", "Tmem184a", "Dntt", "Igll1", "Runx1",
          "Vpreb1", "Rag1", "Rag2", "Xrcc6", "Cd34", "Cd3g", "Cd3e", "Klrb1c", "Tacr1", "Itgam", "Ly6a", "Ly6g", "Ptprc")

meta = readRDS(args$meta)$meta
counts = readRDS(args$expr)$counts
rownames(counts) = idmap$gene(rownames(counts), to="mgi_symbol", dset="mmusculus_gene_ensembl")
counts = counts[rowSums(counts) > 0 & !is.na(rownames(counts)),]
narray::intersect(meta$sample, counts, along=2)

eset = DESeq2::DESeqDataSetFromMatrix(counts, meta, ~1) %>%
    DESeq2::estimateSizeFactors()
reads = DESeq2::counts(eset, normalized=TRUE)
expr = DESeq2::varianceStabilizingTransformation(eset)

saveRDS(list(meta=meta, reads=reads, vst=expr, genes=genes), file=args$outfile)

exp1 = t(assay(expr[genes,]))
meta2 = meta %>% filter(type %in% c("Myeloid", "B-like") | sample == "449s")
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

meta21 = meta2 %>% mutate(cl = clust2$membership) %>% filter(! cl %in% c(1,2))
exp21 = t(assay(expr[genes, meta21$sample]))
pr21 = prcomp(exp21)
clust21 = igraph::cluster_louvain(scran::buildSNNGraph(t(exp21), k=5))
umap21 = uwot::umap(exp21)
colnames(umap21) = c("umap1", "umap2")
umap21 = cbind(meta21, umap21, clust=factor(clust21$membership))

diffg = c("Cd34", "Ptprc", "Ly6a", "Kit")
meta3 = meta2
exp3 = t(assay(expr[diffg, meta3$sample]))
pr3 = prcomp(exp3)

g4 = c("Ets1", "Erg", "Kit")
meta4 = meta #%>% filter(type %in% c("Myeloid", "B-like"))
exp4 = t(assay(expr[g4, meta4$sample]))
pr4 = prcomp(exp4)

pdf(args$plotfile, 10, 8)
plt$pca(pr, aes(x=PC1, y=PC2), annot=meta, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(color=type, size=aneuploidy)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

ggplot(umap, aes(x=umap1, y=umap2)) +
    geom_point(aes(color=type, size=aneuploidy, shape=clust)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()

cbind(umap, t(reads[genes,])) %>% # todo: better to melt reads, then join
    tidyr::gather("gene", "reads", -(sample:clust)) %>%
    ggplot(aes(x=clust, y=log10(reads+1), group=clust)) +
        geom_point(aes(color=type, size=aneuploidy), alpha=0.3) +
        facet_wrap(~ gene, scales="free_y")

plt$pca(pr2, aes(x=PC1, y=PC2), annot=meta2, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(size=aneuploidy, color=type)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

ggplot(umap2, aes(x=umap1, y=umap2)) +
    geom_point(aes(size=aneuploidy, shape=clust, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()

plt$pca(pr21, aes(x=PC1, y=PC2), annot=meta21, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(size=aneuploidy, color=type)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

#ggplot(umap21, aes(x=umap1, y=umap2)) +
#    geom_point(aes(size=aneuploidy, shape=clust, color=type)) +
#    ggrepel::geom_text_repel(aes(label=sample), size=2) +
#    theme_classic()

cbind(umap2, t(reads[genes,meta2$sample])) %>%
    tidyr::gather("gene", "reads", -(sample:clust)) %>%
    ggplot(aes(x=clust, y=log10(reads+1), group=clust)) +
        geom_point(aes(size=aneuploidy, color=type), alpha=0.3) +
        facet_wrap(~ gene, scales="free_y")

plt$pca(pr3, aes(x=PC1, y=PC2), annot=meta3, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(size=aneuploidy, color=type)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()

ggplot(cbind(meta, exp1), aes(x=Ptprc, y=Cd34)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Kit, y=Ly6a)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Ly6a, y=Cd34)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()

plt$pca(pr4, aes(x=PC1, y=PC2), annot=meta4, biplot=TRUE, bi_color="black", bi_size=3) +
    geom_point(aes(size=aneuploidy, color=type)) +
    geom_text(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Ets1, y=Erg)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Kit, y=Erg)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Ly6a, y=Erg)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Ly6a, y=aneuploidy)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    geom_smooth(aes(color=type), method="lm", se=FALSE) +
    theme_classic()
ggplot(cbind(meta, exp1), aes(x=Kit, y=aneuploidy)) +
    geom_point(aes(size=aneuploidy, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    geom_smooth(aes(color=type), method="lm", se=FALSE) +
    theme_classic()
dev.off()
