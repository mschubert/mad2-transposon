library(dplyr)
library(ggplot2)
sys = import('sys')
gset = import('genesets')
util = import('./util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('e', 'eset', 'gene expression rds', 'eset_Mad2PB.rds'),
    opt('d', 'diff_expr', 'rds', 'de_Mad2PB.rds'),
    opt('s', 'setfile', 'rds', '../data/genesets/mouse/MSigDB_Hallmark_2020.rds'),
    opt('o', 'outfile', 'results rds', 'de_Mad2PB/MSigDB_Hallmark_2020.rds'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB/MSigDB_Hallmark_2020.pdf')
)

hl = yaml::read_yaml(args$config)$highlight_de
eset = readRDS(args$eset)$eset
res = readRDS(args$diff_expr)

sets = readRDS(args$setfile) %>%
    gset$filter(min=5, valid=rownames(eset))

pdf(args$plotfile)
for (rname in names(res))
    print(util$plot_gset(res[[rname]], sets) + ggtitle(rname))
dev.off()

saveRDS(res, file=args$outfile)
