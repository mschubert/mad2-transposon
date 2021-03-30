library(dplyr)
library(DESeq2)
sys = import('sys')
util = import('./util')

args = sys$cmd$parse(
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('e', 'eset', 'mad2pb eset', 'eset_Mad2PB.rds'),
    opt('c', 'cis', 'poisson result', '../cis_analysis/poisson.rds'),
    opt('o', 'outfile', 'results rds', 'eset_Mad2PB+multiCIS.rds'),
    opt('p', 'plotfile', 'pdf', 'eset_Mad2PB+multiCIS.pdf')
)

inc_ins = yaml::read_yaml(args$config)$de_ins

cis = readRDS(args$cis)$samples %>%
    filter(external_gene_name %in% inc_ins) %>%
    select(sample, ins=external_gene_name) %>%
    mutate(present = 1) %>%
    narray::construct(present ~ sample + ins, fill=0, data=.)

dset = readRDS(args$eset)
meta = SummarizedExperiment::colData(dset$eset) %>%
    as.data.frame() %>%
    mutate(aneup0.2 = pmin(aneuploidy, 0.2))
keep = meta$sample %in% rownames(cis)
cis = cis[meta$sample[keep],]
idx = cbind(meta[keep,], cis)

vs = dset$vs[,keep]
eset = dset$eset[,keep]
colData(eset) = DataFrame(idx)

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

saveRDS(list(eset=eset, vs=vs, pca=pca, cis=cis, inc_ins=inc_ins), file=args$outfile)
