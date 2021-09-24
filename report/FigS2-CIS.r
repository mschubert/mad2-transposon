library(dplyr)
library(ggplot2)
library(ggraph)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')
bnet = import('../cis_analysis/bionet')

cis_row = function(bn, assocs) {
    p1 = plt$volcano(assocs, size="n_smp", label_top=50, max.overlaps=50, x_label_bias=0.5)
    p2 = ggraph(bn$cis_net) +
        geom_node_point(aes(size=n_smp)) +
        geom_node_text(aes(label=name), repel=TRUE) +
        geom_edge_link(alpha=0.2) +
        theme_void()
    p1 + p2 + plot_layout(widths=c(1,2))
}

aneup_row = function(bn, ext) {
    p1 = plt$volcano(ext, label_top=30)
    net_with_stats = bn$ext_nets$aneuploidy %>%
        left_join(ext %>% select(name=external_gene_name, p.value, statistic))
    p2 = ggraph(net_with_stats) +
        geom_node_point(aes(size=n_smp, fill=statistic, stroke=p.value<0.05,
                                     color=p.value<0.05), shape=21) +
        geom_node_text(aes(label=name), repel=TRUE) +
        geom_edge_link(alpha=0.2) +
        theme_void() +
        scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
        scale_color_manual(name="signif", labels=c("n.s.", "p<0.05"),
                           values=c("white", "black"))
    p1 + p2 + plot_layout(widths=c(1,2))
}

sys$run({
    args = sys$cmd$parse(
#        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
#        opt('p', 'poisson', 'rds', '../cis_analysis/poisson.rds'),
#        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
#        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
#        opt('r', 'rna_ins', 'txt', '../data/rnaseq_imfusion/insertions.txt'),
        opt('p', 'plotfile', 'pdf', 'FigS2-CIS.pdf')
    )

    bn = readRDS("../cis_analysis/bionet_omnipath.rds")
    cis_nodes = igraph::V(bn$cis_net)$name
#    an_nodes = igraph::V(bn$ext_nets$aneuploidy)$name

    assocs = readRDS("../cis_analysis/poisson.rds")$result %>%
        filter(! external_gene_name %in% c("Sfi1", "Drg1")) %>% # excl genome region
        mutate(circle = external_gene_name %in% cis_nodes)
    ext = readRDS("../cis_analysis/ext_gene.rds")$aneuploidy %>%
        mutate(circle = external_gene_name %in% cis_nodes)

    asm = (cis_row(bn, assocs) / aneup_row(bn, ext)) +
        plot_annotation(tag_levels='a') +
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 10, 6)
    print(asm)
    dev.off()
})
