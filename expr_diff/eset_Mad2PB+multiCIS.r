library(dplyr)
library(DESeq2)
io = import('io')
sys = import('sys')
util = import('./util')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'eset', 'mad2pb eset', 'eset_Mad2PB.RData'),
    opt('c', 'cis', 'poisson result', '../cis_analysis/poisson.RData'),
    opt('o', 'outfile', 'results RData', 'eset_Mad2PB+multiCIS.RData'),
    opt('p', 'plotfile', 'pdf', 'eset_Mad2PB+multiCIS.pdf'))

inc_ins = c("Ets1", "Erg", "Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1",
            "Ikzf1", "Trp53", "Rapgef6", "Pten", "Foxn3", "Crebbp",
            "Rpl5", "Cbx5", "Sp3", "Xrcc6", "Notch1")

cis = io$load(args$cis)$samples %>%
    filter(external_gene_name %in% inc_ins) %>%
    select(sample, ins=external_gene_name) %>%
    mutate(present = 1) %>%
    narray::construct(present ~ sample + ins, fill=0, data=.)

dset = io$load(args$eset)
meta = SummarizedExperiment::colData(dset$eset) %>%
    as.data.frame() %>%
    mutate(aneup0.2 = pmin(aneuploidy, 0.2))
keep = meta$sample %in% rownames(cis)
meta = cbind(meta[keep,], cis[meta$sample[keep],])

vs = dset$vs[,keep]
eset = dset$eset[,keep]
eset@colData = DataFrame(meta)

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

save(eset, vs, pca, file=args$outfile)
