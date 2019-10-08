library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
util = import('../expr_diff/util')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('a', 'diff_aneup', 'RData', '../expr_diff/de_Mad2PB.RData'),
    opt('p', 'plotfile', 'pdf', 'sets/KEA_2015.pdf'))

stat1 = readRDS(args$diff_expr)
aneup = io$load(args$diff_aneup)$aneuploidy
aneup$gene_name = toupper(aneup$gene_name)
#aneup$gene_name = idmap$orthologue(aneup$gene_name, dset="mmusculus_gene_ensembl",
#                                   from="mgi_symbol", to="hgnc_symbol")

both = list(stat1=stat1$rev48_stat1_over_wt, aneup=aneup) %>%
    bind_rows(.id="dset") %>%
    filter(padj < 0.1) %>%
    transmute(dset = dset,
              gene_name = gene_name,
              stat = log2FoldChange / lfcSE) %>%
    tidyr::pivot_wider(names_from="dset", values_from="stat") %>%
    na.omit()

ggplot(both, aes(x=stat1, y=aneup)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=gene_name))

pdf(args$plotfile)
for (rname in names(res))
    print(util$plot_gset(res[[rname]], sets) + ggtitle(rname))
dev.off()
