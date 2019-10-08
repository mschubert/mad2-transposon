library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
util = import('../expr_diff/util')

merge_one = function(subs) {
    both = list(stat1=stat1[[subs]], aneup=aneup) %>%
        bind_rows(.id="dset") %>%
        transmute(dset = dset,
                  gene_name = label,
                  padj = adj.p,
                  stat = statistic) %>%
        tidyr::pivot_wider(names_from="dset", values_from=c("stat", "padj")) %>%
        na.omit()
}

plot_one = function(merged) {
    merged = merged %>%
        mutate(label = ifelse(padj_stat1 < 0.1 & padj_aneup < 0.1, gene_name, NA)) %>%
        filter(padj_stat1 < 0.1 | padj_aneup < 0.1)
    ggplot(merged, aes(x=stat_stat1, y=stat_aneup)) +
        geom_point(color="grey") +
        ggrepel::geom_text_repel(aes(label=label), size=3) +
        geom_hline(yintercept=0, linetype="dashed", size=2, alpha=0.3) +
        geom_vline(xintercept=0, linetype="dashed", size=2, alpha=0.3)
}

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('a', 'diff_aneup', 'RData', '../expr_diff/de_Mad2PB.RData'),
    opt('h', 'set_human', 'RData', '../data/genesets/human/KEA_2015.RData'),
    opt('m', 'set_mouse', 'RData', '../data/genesets/mouse/KEA_2015.RData'),
    opt('p', 'plotfile', 'pdf', 'sets/KEA_2015.pdf'))

sets_human = io$load(args$set_human) %>% gset$filter(min=5)
sets_mouse = io$load(args$set_mouse) %>% gset$filter(min=5)

stat1 = readRDS(args$diff_expr) %>%
    lapply(util$test_gsets, set=sets_human)
aneup = io$load(args$diff_aneup)$aneuploidy %>%
    util$test_gsets(set=sets_mouse)

pdf(args$plotfile)
plot_one(merge_one("rev24_cgas_over_wt")) + ggtitle("rev24_cgas_over_wt")
plot_one(merge_one("rev24_stat1_over_wt")) + ggtitle("rev24_stat1_over_wt")
plot_one(merge_one("rev48_stat1_over_wt")) + ggtitle("rev48_stat1_over_wt")
dev.off()
