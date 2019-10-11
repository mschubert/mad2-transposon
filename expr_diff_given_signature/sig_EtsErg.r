library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'other analysis input', '../expr_diff/de_Mad2PB+EtsErg.RData'),
    opt('o', 'outfile', 'rds', 'sig_EtsErg.rds'))

dset = io$load(args$infile)$group_Erg.spleen_vs_Ets1.spleen %>%
    mutate(wald = log2FoldChange / lfcSE) %>%
    head(100) %$%
    setNames(wald, toupper(gene_name))

saveRDS(dset, file=args$outfile)
