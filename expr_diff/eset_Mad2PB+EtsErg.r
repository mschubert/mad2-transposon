library(DESeq2)
library(dplyr)
io = import('io')
sys = import('sys')
util = import('./util')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB.RData'),
    opt('o', 'outfile', 'results RData', 'eset_Mad2PB+EtsErg.RData'),
    opt('p', 'plotfile', 'pdf', 'eset_Mad2PB+EtsErg.pdf'))

eset = io$load(args$eset)$eset
idx = colData(eset) %>%
    as.data.frame() %>%
    mutate(aneup0.3 = pmin(aneuploidy, 0.3),
           keep = type %in% c("T-cell", "Other") & !sample %in% c("401t", "403t", "612t", "631s", "477t"),
           group = factor(DESeq2::counts(eset, normalized=TRUE)["Ets1",] > 700))
levels(idx$group) = c("Erg", "Ets1")
idx$group = paste(idx$group, idx$tissue, sep=":")
idx$group[!idx$keep] = NA
idx$group = relevel(factor(idx$group), "Ets1:spleen")
eset@colData = DataFrame(idx %>% select(-keep))
eset = eset[,idx$keep]
idx = idx[idx$keep,]

design(eset) = ~ group
vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))

pdf(args$plotfile)
pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

save(eset, vs, pca, file=args$outfile)
