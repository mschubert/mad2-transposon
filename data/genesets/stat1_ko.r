library(dplyr)
sys = import('sys')
gset = import('genesets')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', '../../expr_stat1/diff_expr.rds'),
    opt('h', 'human', 'save to rds', 'human/stat1_ko.rds'),
    opt('m', 'mouse', 'save to rds', 'mouse/stat1_ko.rds')
)

# BT549 STAT1 KO differential expression
sets = readRDS(args$diff_expr) %>%
    lapply(. %>% filter(log2FoldChange > 0) %>% head(100) %>% pull(gene_name))

# DoRothEA STAT1 binding + inflammatory MSigDB Hallmarks
gets = gset$get_human(c("MSigDB_Hallmark_2020", "DoRothEA"), conf="a")
msets = c("Interferon Gamma Response", "Interferon Alpha Response",
          "Inflammatory Response", "TNF-alpha Signaling via NF-kB",
          "IL-6/JAK/STAT3 Signaling")
msig = gets$MSigDB_Hallmark_2020[msets] %>% unlist() %>% unique()
sets$`Stat1+IFN` = intersect(msig, gets$DoRothEA$`STAT1 (a)`)

mouse = gset$hu2mouse(sets)

saveRDS(sets, file=args$human)
saveRDS(mouse, file=args$mouse)
