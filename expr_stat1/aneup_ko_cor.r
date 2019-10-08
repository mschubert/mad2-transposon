library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')

merge_one = function(subs) {
    both = list(stat1=stat1[[subs]], aneup=aneup) %>%
        bind_rows(.id="dset") %>%
        transmute(dset = dset,
                  gene_name = gene_name,
                  padj = padj,
                  stat = log2FoldChange / lfcSE) %>%
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
    opt('p', 'plotfile', 'pdf', 'aneup_ko_cor.pdf'))

stat1 = readRDS(args$diff_expr)
aneup = io$load(args$diff_aneup)$aneuploidy
aneup$gene_name = toupper(aneup$gene_name)
#aneup$gene_name = idmap$orthologue(aneup$gene_name, dset="mmusculus_gene_ensembl",
#                                   from="mgi_symbol", to="hgnc_symbol")

pdf(args$plotfile, 10, 8)
plot_one(merge_one("rev24_cgas_over_wt")) + ggtitle("rev24_cgas_over_wt")
plot_one(merge_one("rev24_stat1_over_wt")) + ggtitle("rev24_stat1_over_wt")
plot_one(merge_one("rev48_stat1_over_wt")) + ggtitle("rev48_stat1_over_wt")
dev.off()
