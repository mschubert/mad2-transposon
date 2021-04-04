library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
library(ggraph)
library(ggtext)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')

insertion_matrix = function(cis, rna_ins, cis_aneup, aneup, net_genes) {
    cis_aneup = cis_aneup %>%
        filter(p.value < 0.1)
    rna_ins = rna_ins %>%
        select(sample, external_gene_name=gene_name) %>%
        filter(external_gene_name %in% net_genes) %>%
        distinct() %>%
        mutate(rna_ins = 1)

    cis_result = cis$result %>%
        filter(! external_gene_name %in% c("Sfi1", "Drg1", "Eif4enif1", "Wrap53"), # should blacklist those when mapping
               external_gene_name %in% net_genes) %>%
        filter(adj.p < 0.015)
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
        plot_layout(widths=c(6,1), heights=c(1,0.12,5), guides="collect") &
        theme(plot.margin=margin(0.25, 0, 0.25, 2, "mm"))
}

subtype_assocs = function(ext, net_genes) {
    lvl = setNames(c("Myeloid", "T-ALL", "B-like", "Aneuploidy"),
                   c("Myeloid", "T-cell", "B-like", "aneuploidy"))
    types = bind_rows(ext) %>%
        filter(subset %in% names(lvl),
               external_gene_name %in% net_genes) %>%
        arrange(p.value) %>%
        group_by(subset) %>%
        top_n(4, -p.value)
    types$subset = unname(lvl[types$subset])
    types$subset = factor(types$subset, levels=unname(lvl))

    ggplot(types, aes(x=forcats::fct_reorder(external_gene_name, statistic), y=statistic)) +
        geom_hline(yintercept=0, color="grey", linetype="dashed") +
        geom_bar(aes(fill=subset), stat="identity") +
        geom_text(aes(label=external_gene_name, y=statistic/2)) +
        coord_flip() +
        facet_wrap(~ subset, scales="free_y") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.line = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = FALSE) +
        labs(y = "Wald statistic")
}

bionet_combine = function(bionet) {
    aneup_centrality = bionet$ext_nets$aneuploidy %E>%
        mutate(`Edge type` = "Aneuploidy") %N>%
        select(external_gene_name) %>%
        mutate(aneup_hub = centrality_hub()) %>%
        arrange(-aneup_hub)
    cnet = bionet$cis_net %E>%
        mutate(`Edge type` = "CIS") %N>%
        select(external_gene_name, n_smp) %>%
        mutate(hub = centrality_hub()) %>%
        graph_join(aneup_centrality, copy=TRUE) %>%
        mutate(deg = igraph::degree(.),
               hub = ifelse(is.na(hub), 0, hub)) %>%
        filter(deg > 2) %>%
        arrange(-hub) %>%
        mutate(aneup_hub = ifelse(is.na(aneup_hub), 0, aneup_hub))
    ggraph(cnet) +
        geom_edge_link(aes(color=`Edge type`, alpha=`Edge type`, width=`Edge type`)) +
        scale_edge_color_manual(values=c(CIS="black", Aneuploidy="pink")) +
        scale_edge_width_manual(values=c(CIS=3, Aneuploidy=1)) +
        scale_edge_alpha_manual(values=c(CIS=0.05, Aneuploidy=0.8)) +
        geom_node_point(aes(size=hub, fill=aneup_hub), color="black", shape=21) +
        scale_fill_distiller(palette="RdPu", direction=1) +
        geom_node_label(aes(label=external_gene_name, size=hub), repel=TRUE, min.segment.length=0.75,
            label.size=NA, segment.alpha=0.3, fill="#ffffff00", label.padding=unit(0.31, "lines")) +
        scale_size(range = c(2.5,10)) +
        guides(fill = guide_legend(title="Aneuploidy\ncentrality", override.aes=list(size=5)),
               size = guide_legend(title="CIS centrality"))
}

sc_wgs = function(scs) {
    plt$genome$heatmap(scs) +
        guides(fill = guide_legend(title="Copy number")) +
        theme(panel.ontop = FALSE)
}

aneup_het = function(scs) {
    smp_excl_chr = function(aneuHMM, chr=1:19) {
        lapply(aneuHMM, function(a) {
            a$bins = a$bins[seqnames(a$bins) %in% chr]
            a
        })
    }
    smp_chrs = lapply(scs, smp_excl_chr)
    measures = lapply(smp_chrs, . %>% AneuFinder::karyotypeMeasures() %>% `$`(genomewide)) %>%
        bind_rows(.id="sample")

    m = lm(Heterogeneity ~ 0 + Aneuploidy, data=measures)
    mb = broom::tidy(m)
    mlab = sprintf("p=%.2g<br/>R<sup>2</sup>=%.2f", mb$p.value, broom::glance(m)$r.squared)
    ggplot(measures, aes(x=Aneuploidy, y=Heterogeneity)) +
        geom_abline(slope=mb$estimate, color="purple", linetype="dashed", size=1) +
        geom_point(size=5, alpha=0.6) +
        ggrepel::geom_label_repel(aes(label=sample), box.padding=0.4, segment.color=NA,
                                  label.size=NA, fill="#ffffffc0") +
        annotate("richtext", x=0.1, y=0.25, hjust=0.5, vjust=0.5, label=mlab, fill=NA, label.color=NA) +
        coord_cartesian(clip="off") +
        xlim(0, NA) +
        ylim(0, NA)
}

schema = function() {
    schema = grid::rasterGrob(magick::image_read("external/CISschema.svg"))
    splot = ggplot() +
        annotation_custom(schema) +
        theme(axis.line=element_blank())
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
        select(genotype, sample, type, aneuploidy)
    ext = readRDS(args$ext)
    bionet = readRDS(args$bionet)
    net_genes = bionet$cis_net %N>% pull(external_gene_name)

    # insertion processing
    rna_ins = readr::read_tsv(args$rna_ins)
    cis = readRDS(args$poisson)

    # single-cell shallow WGS
    smps = c("401t", "419t", "413s")
    scs = file.path("../data/wgs", paste0(smps, ".rds")) %>%
        lapply(readRDS) %>%
        setNames(smps)

    # create plot objects
    splot = wrap_plots(wrap_elements(schema() + theme(plot.margin=margin(5,-20,10,-20,"mm"))))
    sc_wgs = wrap_plots(wrap_elements(sc_wgs(scs) + theme(plot.margin = margin(3,0,0,0,"mm"))))
    aneup_het = wrap_plots(wrap_elements(aneup_het(scs) + theme(plot.margin = margin(15,5,15,0,"mm"))))
    stype = wrap_plots(wrap_elements(subtype_assocs(ext, net_genes) + theme(plot.margin = margin(5,5,5,5,"mm"))))

    ins_mat = wrap_plots(wrap_elements(insertion_matrix(cis, rna_ins, ext$aneuploidy, aneup, net_genes) +
                                       theme(plot.margin = margin(0,0,0,0,"cm"))))
    bnet = wrap_plots(wrap_elements(bionet_combine(bionet) + theme(plot.margin = margin(10,0,10,-25,"mm"))))

    top = (splot | sc_wgs | aneup_het | stype) + plot_layout(widths=c(4,5.2,2,3))
    bottom =  wrap_plots(ins_mat) + bnet + plot_layout(widths=c(2.1,1))

    asm = (top / bottom) + plot_layout(heights=c(1,2)) +
        plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig2-CIS.pdf", 20, 12)
    print(asm)
    dev.off()
})
