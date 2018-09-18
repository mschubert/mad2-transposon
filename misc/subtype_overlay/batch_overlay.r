library(dplyr)
library(cowplot)
b = import('base')
io = import('io')
st = import('stats')
idmap = import('process/idmap')
sys = import('sys')

args = sys$cmd$parse(
    opt('m', 'mile', 'expr rdata', '../../data/arrayexpress/E-GEOD-13159.RData'),
    opt('t', 'mad2pb', 'rna-seq rdata', '../../data/rnaseq/assemble.RData'),
    opt('a', 'aneup', 'mile aneup rdata', 'aneup_scores_mile.RData'),
    opt('p', 'plotfile', 'pdf', 'batch_overlay.pdf'))

plot_overlay = function(meta) {
    ggplot(meta, aes(x=x, y=y, shape=factor(covar))) +
        geom_point(aes(fill=type1, color=src, size=src), alpha=0.8) +
        scale_fill_discrete(na.value="#ffeeff00") +
        scale_shape_manual(values=c(21,22,23,24)) +
        scale_color_manual(values=c("black", "#ffeeff00")) +
        scale_size_manual(values=c(3,2)) +
        shadowtext::geom_shadowtext(aes(label=label), size=3, color="black", bg.color="#ffffff55") +
        guides(fill = guide_legend(override.aes=list(shape=21, color="#ffeeff00", size=4)))
}

plot_aneup = function(aneup) {
    ord = aneup %>%
        group_by(type) %>%
        summarize(med = median(aneuploidy)) %>%
        arrange(-med) %>%
        pull(type)
    aneup$type = factor(aneup$type, levels=ord)
    ggplot(aneup, aes(x=type, y=aneuploidy)) +
        geom_violin(color="#00888855", fill="#00000033", alpha=0.2, draw_quantiles=0.5) +
        ggbeeswarm::geom_quasirandom(alpha=0.1) +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=6))
}

mile = io$load(args$mile)
meta_mile = Biobase::pData(mile) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("file") %>%
    transmute(file=file, type1=FactorValue..LEUKEMIA.CLASS., src="mile")

mad2pb = io$load(args$mad2pb)
rownames(mad2pb$expr) = idmap$orthologue(rownames(mad2pb$expr),
    to="hsapiens_homolog_ensembl_gene", dset="mmusculus_gene_ensembl")
mad2pb$expr = mad2pb$expr[!is.na(rownames(mad2pb$expr)),]
meta_mad2pb = mad2pb$idx %>%
    transmute(sample=sample, type2=type, src="mad2pb")

expr = narray::stack(Biobase::exprs(mile), mad2pb$expr, along=2) %>%
    na.omit()
meta = bind_rows(meta_mile, meta_mad2pb) %>%
    mutate(covar = ifelse(is.na(type2), type1, type2),
           label = covar,
           label = ifelse(duplicated(label), NA, label),
           covar = ifelse(grepl("T-cell|T-ALL", covar), "Tcell", covar),
           covar = ifelse(grepl("Myeloid|AML", covar), "Myeloid", covar))
mod = narray::mask(meta$covar, along=2) + 0
corrected = sva::ComBat(dat=expr, batch=factor(meta$src), mod=mod)

tsne = Rtsne::Rtsne(t(corrected), perplexity=50)
pca = prcomp(t(corrected), scale=FALSE)
meta$x = pca$x[,1]
meta$y = pca$x[,2]
meta$covar = ifelse(meta$covar %in% c("Tcell", "Myeloid", "Other"), meta$covar, "Other")
meta$type1 = factor(meta$type1) %>%
    relevel("T-ALL") %>%
    relevel("ALL with t(12;21)") %>%
    relevel("ALL with t(1;19)") %>%
    relevel("ALL with hyperdiploid karyotype")

pdf(args$plotfile, 12, 8)
plot_aneup(io$load(args$aneup)$aneup)

plot_overlay(meta) +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
         title = "Overlay MILE + Mad2PB (PCA)")

meta$x = tsne$Y[,1]
meta$y = tsne$Y[,2]

plot_overlay(meta) +
    labs(x = "t-SNE 1",
         y = "t-SNE 2",
         title = "Overlay MILE + Mad2PB (t-SNE)")

dev.off()
