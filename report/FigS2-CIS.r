library(dplyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')
seq = import('seq')
bnet = import('../cis_analysis/bionet')
tx = seq$coords$transcript(dset="mmusculus_gene_ensembl", assembly="GRCm38")

ins_distr = function(gene, rna_ins, dna_ins) {
    txg = tx %>% filter(external_gene_name == gene)
    rng = txg %>% group_by(ensembl_transcript_id) %>%
        summarize(start=min(transcript_start), end=max(transcript_end))

    exon = ensembldb::exonsBy(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                              filter=AnnotationFilter::GeneIdFilter(txg$ensembl_gene_id[1])) %>%
        as.data.frame() %>% as_tibble() %>%
        select(ensembl_transcript_id=group_name, exon_id, start, end) %>%
        mutate(type="transcript")
    strn = c("1"="+","-1"="-")[as.character(txg$strand[1])]

    ri = rna_ins %>% filter(gene_name == gene, support >= 2) %>%
        mutate(strand=c("1"="+","-1"="-")[as.character(strand)])
    di = dna_ins %>% filter(chr == txg$chromosome_name[1],
                            position >= min(rng$start),
                            position <= max(rng$end)) %>%
        group_by(sample) %>%
            slice_max("reads", n=1) %>%
        ungroup()
    both = bind_rows(list(RNA=ri, DNA=di), .id="type") %>%
        mutate(type = factor(type, levels=c("DNA", "RNA")),
               sense = factor(strand == strn, levels=c("TRUE", "FALSE")))
    levels(both$sense) = c("sense", "antisense")

    p1 = ggplot(both, aes(y=1, x=position)) +
        ggbeeswarm::geom_quasirandom(aes(shape=factor(strand), color=type, alpha=sense), size=7, groupOnX=FALSE) +
        geom_text(aes(label=sample), size=3, position=ggbeeswarm::position_quasirandom(groupOnX=FALSE)) +
        scale_shape_manual(values = c("+"="\u2b95", "-"="\u2b05"), drop=FALSE, guide="none") +
        scale_alpha_manual(values = c("sense"=0.9, "antisense"=0.4), drop=FALSE, name="Orientation") +
        scale_color_discrete(drop=FALSE, name="Evidence") +
        ggtitle(gene) +
        theme(axis.line.y = element_blank(),
              axis.title = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=10)) +
        labs(subtitle = sprintf("chr%s:%.1f-%.1f Mb (%s)", txg$chromosome_name[1],
                                min(rng$start)/1e6, max(rng$end)/1e6, strn))

    p2 = ggplot(txg, aes(y=ensembl_transcript_id, yend=ensembl_transcript_id)) +
        geom_segment(data=rng, aes(x=start, xend=end), size=1, color="#34343456") +
        geom_segment(data=exon, aes(x=start, xend=end), size=5, color="#343434d0") +
        theme_void()

    p1 / p2 + plot_layout(heights=c(1,1))
}

all_ins = function(rna_ins, dna_ins) {
    ps = lapply(c("Ets1", "Erg", "Trp53", "Stat1"), ins_distr, rna_ins=rna_ins, dna_ins=dna_ins)
    wrap_plots(ps) + plot_layout(guides="collect") & theme(plot.margin=unit(c(5,5,5,5),"mm"))
}

cis_row = function(bn, assocs, common) {
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
    wrap_elements(p2) + (p3 + p4 + plot_layout(tag_level="new")) + plot_layout(widths=c(2,1))
}

aneup_row = function(bn, ext, common) {
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
    wrap_elements(p2) + (p3 + p4 + plot_layout(tag_level="new")) + plot_layout(widths=c(2,1))
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
        opt('d', 'dna_ins', 'rds', '../cis_analysis/analysis_set.rds'),
        opt('r', 'rna_ins', 'txt', '../data/rnaseq_imfusion/insertions.txt'),
        opt('s', 'poisson', 'rds', '../cis_analysis/poisson.rds'),
        opt('p', 'plotfile', 'pdf', 'FigS2-CIS.pdf')
    )

    rna_ins = readr::read_tsv(args$rna_ins)
    dna_ins = readRDS(args$dna_ins)

    bn = readRDS(args$bionet)
    bn$cis_net = bn$cis_net %>% mutate(hub = centrality_hub())
    common = intersect(igraph::V(bn$cis_net)$name,
                       igraph::V(bn$ext_nets$aneuploidy)$name)

    assocs = readRDS(args$poisson)$result %>%
        filter(! external_gene_name %in% c("Sfi1", "Drg1")) %>% # excl genome region
        mutate(circle = external_gene_name %in% common)
    ext = readRDS(args$ext)$aneuploidy %>%
        mutate(circle = external_gene_name %in% common)

    volc_cis = plt$volcano(assocs, size="n_smp", label_top=30, max.overlaps=50, x_label_bias=0.5) +
        labs(x="CIS enrichment", y="Adjusted p-value (FDR)")
    volc_aneup = plt$volcano(ext, label_top=30) + labs(x="Aneuploidy difference if inserted", y="P-value")
    cis_row = cis_row(bn, assocs, common)
    aneup_row = aneup_row(bn, ext, common)

    asm = ((volc_cis / volc_aneup) | (cis_row / aneup_row)) +
        plot_layout(widths=c(1,3)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=22, face="bold"),
              axis.title = element_text(size=14),
              axis.text = element_text(size=12))

#    all_ins(rna_ins, dna_ins)

    pdf(args$plotfile, 20, 14)
    print(asm)
    dev.off()
})
