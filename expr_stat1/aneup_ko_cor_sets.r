library(dplyr)
library(ggplot2)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
util = import('../expr_diff/util')

colors = setNames(c("#a50f15", "#006d2c", "#045a8d", "#cccccc", "#8a8a8a"),
                  c("compensated", "hyper-dereg", "inconsistent", "no change", "only 1 dset"))

merge_one = function(subs, wald=1.5) {
    both = list(stat1=stat1[[subs]], aneup=aneup) %>%
        bind_rows(.id="dset") %>%
        transmute(dset = dset,
                  gene_name = label,
                  padj = adj.p,
                  stat = statistic) %>%
        tidyr::pivot_wider(names_from="dset", values_from=c("stat", "padj")) %>%
        mutate(type = case_when(
                is.na(stat_stat1) | is.na(stat_aneup) ~ "only 1 dset",
                abs(stat_stat1) < wald | abs(stat_aneup) < wald ~ "no change",
                stat_stat1 < 0 & stat_aneup < 0 ~ "compensated",
                stat_stat1 > 0 & stat_aneup > 0 ~ "hyper-dereg",
                TRUE ~ "inconsistent"
            ),
            type = factor(type, levels=names(colors))) %>%
        na.omit()
}

plot_one = function(merged) {
    merged = merged %>%
        mutate(top10 = rank(rank(stat_stat1) + rank(stat_aneup)) <= 10,
               bottom10 = rank(-rank(stat_stat1) - rank(stat_aneup)) <= 10,
               inc10 = rank(-abs(rank(stat_stat1) - rank(stat_aneup))) <= 10,
               label = ifelse((top10 | bottom10 | inc10) & padj_stat1 < 0.1 &
                              padj_aneup < 0.1, gene_name, NA))
    ggplot(merged, aes(x=stat_stat1, y=stat_aneup)) +
        geom_point(aes(color=type)) +
        scale_color_manual(values=colors) +
        ggrepel::geom_label_repel(aes(label=label, color=type), size=3,
                                  na.rm=TRUE, segment.alpha=0.3, fill="#ffffffc0",
                                  label.padding=0.1, max.iter=5e4, min.segment.length=0) +
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

pdf(args$plotfile, 16, 14)
plot_one(merge_one("rev24_cgas_over_wt")) + ggtitle("rev24_cgas_over_wt")
plot_one(merge_one("rev24_stat1_over_wt")) + ggtitle("rev24_stat1_over_wt")
plot_one(merge_one("rev48_stat1_over_wt")) + ggtitle("rev48_stat1_over_wt")
dev.off()
