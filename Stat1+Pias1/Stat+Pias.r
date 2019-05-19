io = import('io')
st = import('stats')
plt = import('plot')
enr = import('tools/enrichr')
idmap = import('process/idmap')

# response signatures: Stat1 ChIP (enrichr), Ifn response
expr = io$load("../expr_diff/eset_Mad2PB.RData")
genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1", "Trp53", "Wrap53")
#mile = io$load("../expr_diff/eset_MILE.RData")$expr
#rownames(mile) = idmap$orthologue(rownames(mile), from="hgnc_symbol", to="mgi_symbol")
ee = io$load("../expr_diff/eset_Mad2PB.RData")$vs
stat1 = ee["Stat1",]
rest = ee[-which(rownames(ee) == "Stat1"),]
cc = data.frame(
    gene = rownames(rest),
    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
) %>% arrange(-cor)

encc = io$load("../data/genesets/mouse/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.RData")
chea = io$load("../data/genesets/mouse/ChEA_2016.RData")
go = io$load("../data/genesets/mouse/GO_Biological_Process_2018.RData")
mile = io$load("../data/genesets/mouse/MILE_regulons.RData")
sets = c(go["regulation of complement activation (GO:0030449)"],
         chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
         chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
         encc["STAT3_CHEA"],
         encc["STAT3_ENCODE"],
         STAT1_complement = list(intersect(
            go[["regulation of complement activation (GO:0030449)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
#         STAT1_mile = list(intersect(
#            mile[["STAT1_up"]],
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_cor = list(intersect(
            head(cc$gene, 1000),
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
)
# + MHC, IFN(g), apop(?) -> do GO enrichment, then choose (or: clust?q
gsva = GSVA::gsva(expr$vs, sets)

scores = rbind(expr$vs[genes,], gsva)

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
p6 = ggplot(both, aes(x=Trp53, y=`STAT1_20625510_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p7 = ggplot(both, aes(x=Wrap53, y=`STAT1_20625510_ChIP-Seq_HELA_Human`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))


p8 = ggplot(both, aes(x=Stat1, y=`STAT1_cor`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))
p9 = ggplot(both, aes(x=Pias1, y=`STAT1_cor`)) +
    geom_point(aes(size=aneuploidy, color=type, shape=ins)) +
    geom_text_repel(aes(label=sample))

pdf("Stat+Pias.pdf", 12, 10)
print(p)
corrplot::corrplot(corr)
print(p1)
print(p2)
print(p3)
other = both %>% mutate(Erg = expr$vs["Erg",])
# lmt = st$lm(`STAT1_20625510_ChIP-Seq_HELA_Human` ~ aneuploidy,
#             subsets=other$type, data=other) %>% tbl_df() %>%
#     broom::tidy()
# print(lmt)
ggplot(other, aes(x=`STAT1_20625510_ChIP-Seq_HELA_Human` > 0, y=aneuploidy)) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom(aes(size=Erg, shape=ins, color=type)) +
    ggrepel::geom_text_repel(aes(label=sample), position=ggbeeswarm::position_quasirandom()) +
    facet_wrap(~type, nrow=1, scales="free_y") +
    labs(title = "Stat1 activity with aneuploidy")
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
dev.off()
