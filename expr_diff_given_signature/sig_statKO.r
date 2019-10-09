library(dplyr)
library(magrittr)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'other analysis input', '../expr_stat1/diff_expr.rds'),
    opt('o', 'outfile', 'rds', 'sig_statKO.rds'))

dset = readRDS(args$infile)$rev24_stat1_over_wt %>%
    mutate(wald = log2FoldChange / lfcSE) %>%
    head(100) %$%
    setNames(wald, gene_name)

saveRDS(dset, file=args$outfile)
