library(dplyr)
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', '../../expr_stat1/diff_expr.rds'),
    opt('h', 'human', 'save to RData', 'human/stat1_ko.rds')#,
#    opt('m', 'mouse', 'save to RData', 'mouse/stat1_ko.rds')
)

sets = readRDS(args$diff_expr) %>%
    lapply(. %>% filter(log2FoldChange > 0) %>% head(100) %>% pull(gene_name))

saveRDS(sets, file=args$human)
