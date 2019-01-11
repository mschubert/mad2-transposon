io = import('io')
plt = import('plot')
enr = import('tools/enrichr')

# response signatures: Stat1 ChIP (enrichr), Ifn response
expr = io$load("../expr_diff/eset_Mad2PB.RData")
genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1")

chea = enr$genes("ChEA_2016")
sets = c("STAT1_20625510_ChIP-Seq_HELA_Human",
         "STAT1_17558387_ChIP-Seq_HELA_Human")#,
gsva = GSVA::gsva(expr$vs, chea[sets])

scores = rbind(expr$vs[genes,], scores)

# expression changes with insertions (in different subtypes)
cis = io$load("../cis_analysis/poisson.RData")$samples %>%
    filter(external_gene_name %in% genes, sample %in% colnames(expr$vs)) %>%
    narray::construct(n_ins ~ sample + external_gene_name, data=., fill=0) %>%
    apply(1, function(x) paste(names(x)[x != 0], collapse="+"))
annot = as.data.frame(SummarizedExperiment::colData(expr$eset))
annot$ins = "none"
annot$ins[match(names(cis), annot$sample)] = cis
annot$ins = relevel(factor(annot$ins), "none")

pca = prcomp(t(scores), scale.=TRUE)
p = plt$pca(pca, aes(x=PC1, y=PC2), annot=annot, biplot=TRUE) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))

corr = st$cor$test(t(scores))
rownames(corr) = colnames(corr) = rownames(scores)

both = cbind(annot, t(scores))
p1 = ggplot(both, aes(x=Stat1, y=Pias1)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p2 = ggplot(both, aes(x=`STAT1_17558387_ChIP-Seq_HELA_Human`, y=`STAT1_20625510_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p3 = ggplot(both, aes(x=Stat1, y=`STAT1_20625510_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p4 = ggplot(both, aes(x=Stat1, y=`STAT1_17558387_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p5 = ggplot(both, aes(x=Mb21d1, y=`STAT1_20625510_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))

pdf("Stat+Pias.pdf", 12, 10)
print(p)
corrplot::corrplot(corr)
print(p1)
print(p2)
print(p3)
other = both %>%
    mutate(Erg = expr$vs["Erg",]) %>%
    filter(type == "Other")
lmt = lm(`STAT1_20625510_ChIP-Seq_HELA_Human` ~ aneuploidy, data=other) %>%
    broom::tidy() %>% filter(term == "aneuploidy") # need more pts for wilcox
ggplot(other, aes(x=`STAT1_20625510_ChIP-Seq_HELA_Human` > 0, y=aneuploidy)) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom(aes(size=Erg)) +
    ggrepel::geom_text_repel(aes(label=sample), position=ggbeeswarm::position_quasirandom()) +
    labs(title = "Other Stat1 activity with aneuploidy",
         subtitle = paste("lm p =", round(lmt$p.value, 4)))
print(p4)
print(p5)
dev.off()
