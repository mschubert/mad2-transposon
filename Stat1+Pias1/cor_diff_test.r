io = import('io')
st = import('stats')
plt = import('plot')
enr = import('tools/enrichr')
util = import('../expr_cor/genenet') #TODO: add genenet for each data set [needs common fmt]

#TODO:
# cor pan + subtypes -> 4x1
# cor diff + pvals (3x3 w/ pancan? -> 4x4)
# for Mad2PB + MILE

# delta cor vs delta expr scatter plot

# delta cor w/ ins -> esp. Ets1 ins w/ B-like @low Stat1

#TODOrunning (non-interaction version)
# diff exp w/ types * important_ins

#TODOne
# Fig 2 inserts

# response signatures: Stat1 ChIP (enrichr), Ifn response
expr = io$load("../expr_diff/eset_Mad2PB.RData")
genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1", "Erg", "Ets1", "Ikzf1", "Sox4",
          "Sp3", "Cbx5", "Rpl5", "Crebbp", "Xrcc6", "Myo5b", "Rapgef6", "Pten", "Foxn3", "Trp53")
chea = enr$genes("ChEA_2016")
sets = c("STAT1_20625510_ChIP-Seq_HELA_Human",
         "STAT1_17558387_ChIP-Seq_HELA_Human")
gsva = GSVA::gsva(expr$vs, chea[sets])
ifng = io$load("../expr_sets/gsva_mad2pb/CH.HALLMARK.RData")
ifng = ifng["HALLMARK_INTERFERON_GAMMA_RESPONSE",colnames(gsva),drop=FALSE]
stopifnot(colnames(expr$vs) == colnames(gsva))
scores = rbind(expr$vs[genes,], gsva, ifng[,colnames(gsva),drop=FALSE])

# expression changes with insertions (in different subtypes)
cis = io$load("../cis_analysis/poisson.RData")$samples %>%
    filter(external_gene_name %in% genes, sample %in% colnames(expr$vs)) %>%
    narray::construct(n_ins ~ sample + external_gene_name, data=., fill=0) %>%
    apply(1, function(x) paste(names(x)[x != 0], collapse="+"))
annot = as.data.frame(SummarizedExperiment::colData(expr$eset))
annot$ins = "none"
annot$ins[match(names(cis), annot$sample)] = cis
annot$ins = relevel(factor(annot$ins), "none")

dset = list(
    `pan` = t(scores),
    `B-like` = t(scores[,annot$type=="Other"]),
    `T-ALL` = t(scores[,annot$type=="T-cell"]),
    `Myeloid` = t(scores[,annot$type=="Myeloid"])
)
corplots = lapply(dset, util$plot_cor_matrix)

# pcor for subtypes are too unstable
pdf("cor_diff_test.pdf", 12, 10)
corrplot::corrplot(cor(dset$`pan`), main="pan")
util$plot_pcor_net(util$pcor(dset$`pan`))
util$plot_bootstrapped_pcor(dset$`pan`, n=1000, fdr=0.3, show_edge_if=100)
corrplot::corrplot(cor(dset$`B-like`), main="B-like")
corrplot::corrplot(cor(dset$`T-ALL`), main="T-ALL")
corrplot::corrplot(cor(dset$`Myeloid`), main="Myeloid")
dev.off()
