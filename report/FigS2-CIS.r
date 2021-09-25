library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')
bnet = import('../cis_analysis/bionet')

cis_row = function(bn, assocs, common) {
    p1 = plt$volcano(assocs, size="n_smp", label_top=30, max.overlaps=50, x_label_bias=0.5) +
        labs(x="CIS enrichment", y="Adjusted p-value (FDR)")
    p2 = ggraph(bn$cis_net) +
        geom_node_point(aes(size=n_smp, alpha=hub)) +
        geom_node_label(aes(label=name), label.size=NA, fill="#ffffffa0",
                        label.padding = unit(0.12, "lines"), repel=TRUE) +
        geom_edge_link(alpha=0.2) +
        theme_void() +
        labs(size = "Samples\nwith\ninsertions", alpha="Hub\ncentrality")
    stat_df = bn$cis_net %N>% as_tibble()
    top_smp = stat_df %>% filter(rank(-n_smp, ties.method="first") <= 25)
    top_hub = stat_df %>% filter(rank(-hub, ties.method="first") <= 25)
    p3 = ggplot(top_smp, aes(x=forcats::fct_reorder(name, n_smp), y=n_smp)) +
        geom_col(aes(alpha=name %in% common)) +
        coord_flip() +
        scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.3), guide=FALSE) +
        labs(x="", y="Samples with insertions")
    p4 = ggplot(top_hub, aes(x=forcats::fct_reorder(name, hub), y=hub)) +
        geom_col(aes(alpha=name %in% common)) +
        coord_flip() +
        scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.3), guide=FALSE) +
        labs(x="CIS gene", y="Hub centrality")
    p1 + p2 + p3 + p4 + plot_layout(widths=c(1.2,2,0.5,0.5))
}

aneup_row = function(bn, ext, common) {
    p1 = plt$volcano(ext, label_top=30) + labs(x="Aneuploidy difference if inserted", y="P-value")
    net_with_stats = bn$ext_nets$aneuploidy %>%
        left_join(ext %>% select(name=external_gene_name, p.value, statistic))
    p2 = ggraph(net_with_stats) +
        geom_node_point(aes(size=n_smp, fill=statistic, stroke=p.value<0.05,
                                     color=p.value<0.05), shape=21) +
        geom_node_label(aes(label=name), label.size=NA, fill="#ffffffa0",
                        label.padding = unit(0.12, "lines"), repel=TRUE) +
        geom_edge_link(alpha=0.2) +
        theme_void() +
        scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
        scale_color_manual(name="Aneuploidy\nsignificance", labels=c("n.s.", "p<0.05"),
                           values=c("white", "black")) +
        labs(fill = "Aneuploidy\nWald\nstatistic",
             size = "Samples\nwith\ninsertions")
    stat_df = bn$ext_nets$aneuploidy %>% mutate(hub = centrality_hub()) %N>% as_tibble()
    top_smp = stat_df %>% filter(rank(-n_smp, ties.method="first") <= 25)
    top_hub = stat_df %>% filter(rank(-hub, ties.method="first") <= 25)
    p3 = ggplot(top_smp, aes(x=forcats::fct_reorder(name, n_smp), y=n_smp)) +
        geom_col(aes(alpha=name %in% common)) +
        coord_flip() +
        scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.3), guide=FALSE) +
        labs(x="", y="Samples with insertions")
    p4 = ggplot(top_hub, aes(x=forcats::fct_reorder(name, hub), y=hub)) +
        geom_col(aes(alpha=name %in% common)) +
        coord_flip() +
        scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0.3), guide=FALSE) +
        labs(x="CIS gene", y="Hub centrality")
    p1 + p2 + p3 + p4 + plot_layout(widths=c(1.2,2,0.5,0.5))
}

sys$run({
    args = sys$cmd$parse(
#        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('p', 'poisson', 'rds', '../cis_analysis/poisson.rds'),
        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
#        opt('r', 'rna_ins', 'txt', '../data/rnaseq_imfusion/insertions.txt'),
        opt('p', 'plotfile', 'pdf', 'FigS2-CIS.pdf')
    )

    bn = readRDS(args$bionet)
    bn$cis_net = bn$cis_net %>% mutate(hub = centrality_hub())
    common = intersect(igraph::V(bn$cis_net)$name,
                       igraph::V(bn$ext_nets$aneuploidy)$name)

    assocs = readRDS(args$poisson)$result %>%
        filter(! external_gene_name %in% c("Sfi1", "Drg1")) %>% # excl genome region
        mutate(circle = external_gene_name %in% common)
    ext = readRDS(args$ext)$aneuploidy %>%
        mutate(circle = external_gene_name %in% common)

    asm = (cis_row(bn, assocs, common) / aneup_row(bn, ext, common)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 18, 14)
    print(asm)
    dev.off()
})
