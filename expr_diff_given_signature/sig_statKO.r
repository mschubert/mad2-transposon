library(dplyr)
library(magrittr)
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('i', 'infile', 'other analysis input', '../expr_stat1/diff_expr.rds'),
    opt('o', 'outfile', 'rds', 'sig_statKO.rds'))

dset = readRDS(args$infile)$rev24_stat1_over_wt %>%
    mutate(gene_name = #idmap$orthologue(gene_name, to="mgi_symbol"),
                       stringr::str_to_title(gene_name), #FIXME: do actual mapping here
           wald = log2FoldChange / lfcSE) %>%
    head(100) %$%
    setNames(wald, gene_name)

saveRDS(dset, file=args$outfile)
