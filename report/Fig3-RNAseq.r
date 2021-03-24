library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
plt = import('plot')

umap_types = function(em) {
    set.seed(15209)
    vst = as.matrix(SummarizedExperiment::assay(em$vst))
    umap2 = uwot::umap(t(vst), n_components=2)
    dset = em$meta %>%
        mutate(umap1 = umap2[,1],
               umap2 = umap2[,2])

    p1 = ggplot(dset, aes(x=umap1, y=umap2)) +
        geom_point(aes(color=type, size=aneuploidy), alpha=0.8) +
        ggrepel::geom_text_repel(aes(label=sample), size=3, segment.alpha=0.3)

    markers = c("Cd3g", "Ly6g", "Ebf1", "Ighm", "Kit", "Ets1", "Erg", "Stat1")
    reads = em$reads[markers,]
    reads["Ighm",] = reads["Ighm",] / 100
#    reads["Ly6g",] = reads["Ly6g",] * 10
    reads = reads %>%
        reshape2::melt() %>%
        as_tibble() %>%
        dplyr::rename(gene=Var1, sample=Var2, reads=value) %>%
        group_by(gene) %>%
            mutate(reads_scaled = reads / max(reads)) %>%
        ungroup() %>%
        inner_join(dset %>% select(sample, umap1, umap2, aneuploidy))
    p2 = ggplot(reads, aes(x=umap1, y=umap2)) +
        geom_point(aes(fill=reads, size=reads_scaled), color="black", alpha=0.8, shape=21) +
        facet_wrap(~gene, nrow=2) +
        scale_fill_distiller(palette="Spectral", trans="log10") +
        coord_cartesian(clip="off") +
        theme(strip.background = element_rect(fill="#efefef"))
    p1 + p2 + plot_layout(nrow=1, widths=c(1.1,2))
}

aneup_volcano = function(diff_expr) {
    diff_expr = diff_expr %>%
        lapply(. %>% mutate(stat = log2FoldChange / lfcSE))
    sets = gset$get_mouse(c("MSigDB_Hallmark_2020", "DoRothEA"), conf="a")
    names(sets[[2]]) = sub(" (a)", "", names(sets[[2]]), fixed=TRUE)
    ifn = unlist(sets[[1]][c("Interferon Gamma Response", "Interferon Alpha Response",
                             "Inflammatory Response", "TNF-alpha Signaling via NF-kB",
                             "IL-6/JAK/STAT3 Signaling")],
                 use.names=FALSE)
    sets[[2]] = c(sets[[2]], "Stat1+IFN"=list(intersect(sets[[2]]$STAT1, ifn)))

    hmFC = gset$test(diff_expr$aneuploidy, sets[[1]], stat="log2FoldChange")
    goFC = gset$test(diff_expr$aneuploidy, sets[[2]], stat="log2FoldChange")
    lfc = bind_rows(hmFC, goFC) %>% select(label, mean_lfc=estimate)

    hm = gset$test(diff_expr$aneuploidy, sets[[1]])
    go = gset$test(diff_expr$aneuploidy, sets[[2]])
    both = bind_rows(hm, go) %>%
        arrange(adj.p, p.value) %>%
        inner_join(lfc) %>%
        plt$p_effect("adj.p", "mean_lfc", thresh=0.01)
    plt$volcano(both, label_top=20, repel=TRUE, max.overlaps=Inf,
                pos_label_bias=1.5, x_label_bias=0.4, base.size=0.5, text.size=4) +
        xlab("Mean log2 fold change in set") +
        ylab("adjusted p-value (FDR) Wald change)") +
        theme(#axis.line.y = element_blank(),
              panel.grid.major = element_line(color="#efefef", size=0.5))
}

gset_cor = function(meta) {
    meta2 = meta %>%
        filter(!is.na(type))

    p1 = ggplot(meta2, aes(x=`Interferon Gamma Response`, y=aneuploidy)) +
        geom_point(aes(color=type, size=STAT1), alpha=0.5) +
        geom_smooth(aes(color=type), method="lm", se=FALSE) +
        geom_smooth(method="lm", se=FALSE, color="black") +
        plot_layout(tag_level="new")
    d1 = ggplot(meta2, aes(x=type, y=`Interferon Gamma Response`)) +
        geom_boxplot(aes(fill=type)) +
        coord_flip() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none")
    p2 = ggplot(meta2, aes(x=`Myc Targets V1`, y=aneuploidy)) +
        geom_point(aes(color=type, size=STAT1), alpha=0.5) +
        geom_smooth(aes(color=type), method="lm", se=FALSE) +
        geom_smooth(method="lm", se=FALSE, color="black") +
        theme(axis.title.y = element_blank()) +
        plot_layout(tag_level="new")
    d2 = ggplot(meta2, aes(x=type, y=`Myc Targets V1`)) +
        geom_boxplot(aes(fill=type)) +
        coord_flip() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none")
    d3 = ggplot(meta2, aes(x=type, y=aneuploidy)) +
        geom_boxplot(aes(fill=type)) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none") +
        plot_layout(tag_level="new")

    d1 + d2 + plot_spacer() + p1 + p2 + d3 +
        plot_layout(widths=c(5,5,1), heights=c(1,5), guides="collect")
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('c', 'cis', 'rds', '../cis_analysis/poisson.rds'),
        opt('e', 'expr', 'rds', '../expr_diff/de_Mad2PB.rds'),
        opt('h', 'gsva_hm', 'rds', '../data/gsva/mad2pb/MSigDB_Hallmark_2020.rds'),
        opt('d', 'gsva_dorothea', 'rds', '../data/gsva/mad2pb/DoRothEA.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig3-RNAseq.pdf')
    )

    ghm = readRDS(args$gsva_hm)
    gdo = readRDS(args$gsva_dorothea)
    narray::intersect(ghm, gdo, along=2)
    sets = tibble(sample = colnames(ghm),
                  `Myc Targets V1` = ghm["Myc Targets V1",],
                  `Interferon Gamma Response` = ghm["Interferon Gamma Response",],
                  STAT1 = gdo["STAT1 (a)",])

    cis = readRDS(args$cis)$samples # todo: add stat1 ins to lm plot

    aset = readRDS(args$aset)
    meta = aset$meta %>%
        left_join(sets)

    markers = readRDS("../expr_markers/markers.rds")
    diff_expr = readRDS(args$expr)

    umap = umap_types(markers)
    volc2 = aneup_volcano(diff_expr)
    gsc = gset_cor(meta)

    asm = umap / (volc2 + (gsc / plot_spacer() + plot_layout(heights=c(1.2,1))) + plot_layout(widths=c(1,1.5))) +
        plot_layout(heights=c(1,2)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 14)
    print(asm)
    dev.off()
})
