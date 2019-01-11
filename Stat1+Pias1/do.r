io = import('io')
plt = import('plot')
enr = import('tools/enrichr')

# expression of Cgas, Stat1, Pias1, Ifng, Ifngr1 in different subtypes + aneup (PCA biplot for genes?)
# expression changes with insertions (in different subtypes)
expr = io$load("../expr_diff/eset_Mad2PB.RData")
genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1")
cis = io$load("../cis_analysis/poisson.RData")$samples %>%
    filter(external_gene_name %in% genes, sample %in% colnames(expr$vs)) %>%
    narray::construct(n_ins ~ sample + external_gene_name, data=., fill=0) %>%
    apply(1, function(x) paste(names(x)[x != 0], collapse="+"))

annot = as.data.frame(SummarizedExperiment::colData(expr$eset))
annot$ins = "none"
annot$ins[match(names(cis), annot$sample)] = cis
annot$ins = relevel(factor(annot$ins), "none")

pca = prcomp(t(expr$vs[genes,]))
p1 = plt$pca(pca, aes(x=PC1, y=PC2), annot=annot, biplot=TRUE) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample)) +
    coord_fixed()

p2 = plt$pca(pca, aes(x=PC3, y=PC4), annot=annot, biplot=TRUE) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample)) +
    coord_fixed()

p3 = plt$pca(pca, aes(x=PC5, y=PC6), annot=annot, biplot=TRUE) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample)) +
    coord_fixed()

pdf("Stat+Pias.pdf", 10, 10)
print(p1)
print(p2)
print(p3)
dev.off()

# response signatures: Stat1 ChIP (enrichr), Ifn response
sets = enr$genes("ChEA_2016")
