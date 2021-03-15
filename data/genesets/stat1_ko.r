library(dplyr)
sys = import('sys')
gset = import('genesets')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', '../../expr_stat1/diff_expr.rds'),
    opt('h', 'human', 'save to rds', 'human/stat1_ko.rds'),
    opt('m', 'mouse', 'save to rds', 'mouse/stat1_ko.rds')
)

sets = readRDS(args$diff_expr) %>%
    lapply(. %>% filter(log2FoldChange > 0) %>% head(100) %>% pull(gene_name))

mouse = gset$hu2mouse(sets)

saveRDS(sets, file=args$human)
saveRDS(mouse, file=args$mouse)
