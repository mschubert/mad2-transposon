library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
library(ggraph)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')

insertion_matrix = function(cis, rna_ins, aneup, net_genes) {
    rna_ins = rna_ins %>%
        select(sample, external_gene_name=gene_name) %>%
        filter(external_gene_name %in% net_genes) %>%
        distinct() %>%
        mutate(rna_ins = 1)

    cis_result = cis$result %>%
        filter(! external_gene_name %in% c("Sfi1", "Drg1", "Eif4enif1"), # should blacklist those when mapping
               external_gene_name %in% net_genes) %>%
        filter(adj.p < 0.01)
    cis_samples = cis$samples %>%
        filter(external_gene_name %in% cis_result$external_gene_name) %>%
        tidyr::complete(sample, external_gene_name, fill=list(n_ins=0, reads=0)) %>%
        group_by(sample) %>%
        mutate(total_ins = sum(n_ins),
               total_reads = sum(reads),
               gene_read_frac = reads / total_reads) %>%
        ungroup() %>%
        left_join(rna_ins) %>%
        mutate(has_ins = ifelse(rna_ins | n_ins != 0, TRUE, NA),
            rna_ins = ifelse(is.na(rna_ins), 0, 1),
            ins_type = case_when(
                n_ins > 0 & rna_ins > 0 ~ "both",
                n_ins > 0 & rna_ins == 0 ~ "DNA",
                n_ins == 0 & rna_ins > 0 ~ "RNA"
            ),
            ins_type = factor(ins_type, levels=c("DNA", "RNA", "both"))
        ) %>%
        inner_join(aneup) %>%
        inner_join(ext$aneuploidy %>% select(external_gene_name, aneup_stat=statistic)) %>%
        mutate(sample = forcats::fct_reorder(sample, aneuploidy),
               external_gene_name = forcats::fct_reorder(external_gene_name, aneup_stat))
    genelvl = levels(cis_samples$external_gene_name)
    smplvl = levels(cis_samples$sample)
    cis_stats = cis_result %>%
        mutate(external_gene_name = factor(external_gene_name, levels=genelvl)) %>%
        filter(!is.na(external_gene_name))

    p1 = ggplot(cis_samples, aes(x=sample, y=external_gene_name)) +
        geom_tile(aes(fill=ins_type, alpha=gene_read_frac, color=has_ins), size=0.2) +
        scale_fill_manual(values=c("maroon4", "navy", "springgreen4"), na.translate=FALSE) +
        scale_color_manual(values="#565656ff") +
        guides(color = FALSE,
               fill = guide_legend(title="Insert type"),
               alpha = guide_legend(title="Read fraction")) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
        labs(x = "Sample",
             y = "Transposon-inserted gene") +
        coord_cartesian(expand=FALSE)

    p11 = mutate(aneup, sample=factor(sample, smplvl)) %>%
        filter(!is.na(sample)) %>%
        ggplot(aes(x=sample, y=aneuploidy, fill=type)) +
            geom_bar(stat="identity") +
            geom_hline(yintercept=0.2, color="grey", linetype="dotted") +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank()) +
            guides(fill=guide_legend(title="Cancer type")) +
            labs(y = "Aneuploidy") +
            theme(plot.tag.position = c(0, 1.2)) +
            coord_cartesian(expand=FALSE)

    p12 = mutate(aneup, sample=factor(sample, smplvl)) %>%
        filter(!is.na(sample)) %>%
        ggplot(aes(x=sample, y=1)) +
            geom_tile(aes(fill=genotype), color="black", size=0.2) +
            theme(axis.text = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank()) +
            guides(fill=guide_legend(title="Genotype")) +
            scale_fill_manual(values=c("Mad2 PB Mx1-Cre"="darkorchid", "PB Mx1-Cre"="white")) +
            coord_cartesian(expand=FALSE)

    p13 = ggplot(cis_stats, aes(x=external_gene_name, y=-log10(adj.p))) +
        geom_bar(stat="identity") +
        scale_x_discrete(position="top") +
        coord_flip(expand=FALSE) +
        geom_hline(yintercept=-log10(0.01), color="grey", linetype="dotted") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank()) +
        labs(y = "-log FDR")

    p11 + plot_spacer() + (p12 + plot_layout(tag_level="new")) + plot_spacer() +
        p1 + (p13 + plot_layout(tag_level="new")) +
        plot_layout(widths=c(6,1), heights=c(1,0.1,5), guides="collect") &
        theme(plot.margin=margin(0.25, 0, 0.25, 2, "mm"))
}

subtype_assocs = function(ext, net_genes) {
    lvl = setNames(c("Myeloid", "T-ALL", "B-like", "Aneuploidy"),
                   c("Myeloid", "T-cell", "Other", "aneuploidy"))
    types = bind_rows(ext) %>%
        filter(subset %in% names(lvl),
               external_gene_name %in% net_genes) %>%
        arrange(p.value) %>%
        group_by(subset) %>%
        top_n(4, -p.value)
    types$subset = unname(lvl[types$subset])
    types$subset = factor(types$subset, levels=unname(lvl))

    p = ggplot(types, aes(x=forcats::fct_reorder(external_gene_name, statistic), y=statistic)) +
        geom_hline(yintercept=0, color="grey", linetype="dashed") +
        geom_bar(aes(fill=subset), stat="identity") +
        geom_text(aes(label=external_gene_name, y=statistic/2)) +
        coord_flip() +
        facet_wrap(~ subset, scales="free_y", ncol=1) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.line = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = FALSE) +
        labs(y = "Wald statistic")
    wrap_plots(wrap_elements(p))
}

bionet_combine = function(bionet) {
    aneup_centrality = bionet$ext_nets$aneuploidy %N>%
        select(external_gene_name) %>%
        mutate(aneup_hub = centrality_hub()) %>%
        arrange(-aneup_hub)
    cnet = bionet$cis_net %N>%
        select(external_gene_name, n_smp) %>%
        mutate(deg = igraph::degree(.),
               hub = centrality_hub()) %>%
        filter(deg > 1) %>%
        mutate(deg = igraph::degree(.)) %>%
        filter(deg > 1) %>%
        arrange(-hub) %>%
        left_join(aneup_centrality, copy=TRUE) %>%
        mutate(aneup_hub = ifelse(is.na(aneup_hub), 0, aneup_hub))
    p = ggraph(cnet) +
        geom_edge_link(alpha=0.05, width=3) +
        geom_node_point(aes(size=hub, fill=aneup_hub), color="black", shape=21) +
        scale_fill_distiller(palette="RdPu", direction=1) +
        geom_node_label(aes(label=external_gene_name, size=hub), repel=TRUE,
                        label.size=NA, segment.alpha=0.3, fill="#ffffff00") +
        scale_size(range = c(2.5,12)) +
        guides(fill = guide_legend(title="Aneuploidy\ncentrality", override.aes=list(size=5)),
               size = guide_legend(title="CIS centrality"))
    wrap_plots(p)
}

sc_wgs = function() {
    smps = c("401t", "419t", "413s")
    scs = file.path("../data/wgs", paste0(smps, ".rds")) %>%
        lapply(readRDS) %>%
        setNames(smps)
    plt$genome$heatmap_aneuHMM(scs) + guides(fill = guide_legend(title="Copy number"))
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('p', 'poisson', 'rds', '../cis_analysis/poisson.rds'),
        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
        opt('r', 'rna_ins', 'txt', '../data/rnaseq_imfusion/insertions.txt')
    )

    aneup = readRDS(args$aset)$meta %>%
        select(genotype, sample, type, aneuploidy) %>%
        mutate(type = factor(type))
    levels(aneup$type)[levels(aneup$type) == "Other"] = "B-like"
    ext = readRDS(args$ext)
    bionet = readRDS(args$bionet)
    net_genes = bionet$cis_net %N>% pull(external_gene_name)

    # insertion processing
    rna_ins = readr::read_tsv(args$rna_ins)
    cis = readRDS(args$poisson)

    # create plot objects
    schema = grid::rasterGrob(magick::image_read("external/CISschema.svg"))
    splot = ggplot() + annotation_custom(schema) +
        theme(axis.line=element_blank())
    sc_wgs = sc_wgs()
    ins_mat = insertion_matrix(cis, rna_ins, aneup, net_genes)
    stype = subtype_assocs(ext, net_genes)
    bnet = bionet_combine(bionet)

    asm = ((splot + sc_wgs + plot_layout(widths=c(2,3))) /
        ins_mat /
        (stype + bnet + plot_layout(widths=c(2,7)))) +
            plot_layout(heights=c(1.2,3,2)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig2-CIS.pdf", 14, 19)
    print(asm)
    dev.off()
})
