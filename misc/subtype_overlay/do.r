library(dplyr)
library(cowplot)
b = import('base')
io = import('io')
st = import('stats')
idmap = import('process/idmap')

mile = io$load("../../data/arrayexpress/E-GEOD-13159.RData")
meta_mile = Biobase::pData(mile) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("file") %>%
    transmute(file=file, type1=FactorValue..LEUKEMIA.CLASS., src="mile")

mad2pb = io$load("../../data/rnaseq/assemble.RData")
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
corrected = sva::ComBat(expr, batch=factor(meta$src), mod=mod)

#tsne = Rtsne::Rtsne(t(corrected))
#meta$x = tsne$Y[,1]
#meta$y = tsne$Y[,2]

pca = prcomp(t(expr), scale=FALSE)
meta$x = pca$x[,2]
meta$y = pca$x[,3]

pdf("subtype_overlay.pdf", 12, 8)
ggplot(meta, aes(x=x, y=y)) +
    geom_point(aes(color=type1), alpha=0.5, shape=16, size=2) +
    geom_text(aes(color=type1, label=label), size=3) +
    geom_point(aes(shape=type2), alpha=1, size=2) +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
         title = "Overlay MILE + Mad2PB")
dev.off()
