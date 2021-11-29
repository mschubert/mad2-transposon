library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
library(ggtext)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
plt = import('plot')

umap_types = function(em) {
    set.seed(15202)
    meta = em$meta %>% filter(analysis_set)
    vst = as.matrix(SummarizedExperiment::assay(em$vst[,meta$sample]))
    umap2 = uwot::umap(t(vst), n_components=2)
    dset = meta %>%
        mutate(umap1 = umap2[,1],
               umap2 = umap2[,2])

    p1 = ggplot(dset, aes(x=umap1, y=umap2)) +
        geom_point(aes(color=type, size=aneuploidy), alpha=0.8) +
        ggrepel::geom_text_repel(aes(label=sample), size=3, segment.alpha=0.3) +
        labs(color="Type", size="Aneuploidy")

    markers = c("Cd3g", "Ly6g", "Lyz1", "Ebf1", "Kit", "Ets1", "Erg", "Stat1")
    reads = em$reads[markers,]
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
        labs(fill="Read count", size="Read count\n(scaled)") +
        theme(strip.background = element_rect(fill="#efefef"))
    p1 + p2 + plot_layout(nrow=1, widths=c(1.1,2))
}

aneup_volcano = function(diff_expr) {
    diff_expr = diff_expr %>%
        lapply(. %>% mutate(stat = log2FoldChange / lfcSE))
    sets = gset$get_mouse(c("MSigDB_Hallmark_2020", "DoRothEA"), conf="a")
    hmFC = gset$test_lm(diff_expr$aneuploidy, sets[[1]], stat="log2FoldChange")
    tfFC = gset$test_lm(diff_expr$aneuploidy, sets[[2]], stat="log2FoldChange")
    lfc = bind_rows(hmFC, tfFC) %>% select(label, mean_lfc=estimate)

    hm = gset$test_lm(diff_expr$aneuploidy, sets[[1]])
    go = gset$test_lm(diff_expr$aneuploidy, sets[[2]])
    both = bind_rows(hm, go) %>%
        arrange(adj.p, p.value) %>%
        inner_join(lfc) %>%
        plt$p_effect("adj.p", "mean_lfc", thresh=0.01)
    plt$volcano(both, label_top=20, repel=TRUE, max.overlaps=Inf,
                pos_label_bias=1.5, x_label_bias=0.4, base.size=0.5, text.size=4) +
        xlab("Mean log2 fold change in set") +
        ylab("adjusted p-value (FDR) Wald change") +
        theme(axis.title = element_text(size=14),
              axis.text = element_text(size=12),
              panel.grid.major = element_line(color="#efefef", size=0.5))
}

gset_aneup = function(dset) {
    cols = c("Myeloid"="#f8766d", "T-cell"="#619cff", #"B-like"="#00ba38",
             "Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen1")
    gdf = dset %>%
        mutate(Subtype = ifelse(is.na(Driver), as.character(Type), as.character(Driver))) %>%
        filter(Subtype != "B-like") # rare B-like

    cor_p = . %>% lm(data=gdf) %>% broom::tidy() %>% filter(term == "Aneuploidy0.2")
    m1n = cor_p(`Interferon Gamma Response` ~ Aneuploidy0.2)
    m1c1 = cor_p(`Interferon Gamma Response` ~ Subtype + Aneuploidy0.2)
    m1c2 = cor_p(`Interferon Gamma Response` ~ Subtype + `STAT1 (a)` + Aneuploidy0.2)
    m2n = cor_p(`Myc Targets V1` ~ Aneuploidy0.2)
    m2c1 = cor_p(`Myc Targets V1` ~ Myc_copies + Aneuploidy0.2)
    m2c2 = cor_p(`Myc Targets V1` ~ Myc_copies + Subtype + Aneuploidy0.2)
    m2c3 = cor_p(`Myc Targets V1` ~ Myc_copies + `STAT1 (a)` + Aneuploidy0.2)
    m1 = sprintf("p=%.2g (naive)<br/>%.2g (Type/Driver)<br/>%.2g (Type/Driver + STAT1)",
                 m1n$p.value, m1c1$p.value, m1c2$p.value)
    m2 = sprintf("p=%.2g (naive)<br/>%.2g (Myc copies)<br/>%.2g (Myc copies + Type/Driver)<br/>%.2g (Myc copies + STAT1)",
                 m2n$p.value, m2c1$p.value, m2c2$p.value, m2c3$p.value)

    common = list(
        geom_point(aes(fill=Subtype, size=`STAT1 (a)`, shape=CIS), color="black", alpha=0.5),
        geom_smooth(aes(color=Subtype), method="lm", se=FALSE),
        scale_fill_manual(values=cols, name="Type/Driver"),
        scale_color_manual(values=cols, name="Type/Driver"),
        scale_shape_manual(values=c("Stat1"=25, "Pias1"=24), na.translate=TRUE, na.value=21, guide=FALSE),
        ylab("Aneuploidy"),
        theme(plot.margin = margin(0,0,0,0,"mm")),
        plot_layout(tag_level="new")
    )
    dens_common = list(
        geom_density(size=0.5, alpha=0.1, adjust=1.5),
        scale_color_manual(values=cols, guide=FALSE),
        scale_fill_manual(values=cols, guide=FALSE),
        theme(plot.margin = margin(0,0,0,0,"mm")),
        theme_void()
    )
    p1 = ggplot(gdf, aes(x=`Interferon Gamma Response`, y=Aneuploidy0.2)) + common +
        annotate("richtext", x=0, y=0.104, label=m1, hjust=0.7, vjust=0.5, angle=-34, size=4,
                 label.size=NA, fill="#ffffff90") +
        geom_smooth(method="lm", se=FALSE, color="black")
    p2 = ggplot(gdf, aes(x=`Myc Targets V1`, y=Aneuploidy0.2)) + common +
        annotate("richtext", x=-0.2, y=0.108, label=m2, hjust=0.4, vjust=0.58, angle=24, size=4,
                 label.size=NA, fill="#ffffffc0") +
        geom_smooth(method="lm", se=FALSE, color="black") +
        theme(axis.title.y = element_blank(),
              plot.margin = margin(0,0,0,0,"mm"))
    d1 = ggplot(gdf, aes(x=`Interferon Gamma Response`, color=Subtype, fill=Subtype)) + dens_common
    d2 = ggplot(gdf, aes(x=`Myc Targets V1`, color=Subtype, fill=Subtype)) + dens_common
    d3 = ggplot(gdf, aes(x=Aneuploidy0.2, color=Subtype, fill=Subtype)) + dens_common +
        coord_flip() + ylab("Aneuploidy") + plot_layout(tag_level="new")

    d1 + d2 + plot_spacer() + p1 + p2 + d3 +
        plot_layout(widths=c(5,5,1), heights=c(1,5), guides="collect")
}

set_tissue = function(dset) {
    tsets = dset %>% filter(!is.na(Type))
    subsets = dset %>% filter(!is.na(Driver))
    keep = c("Interferon Gamma Response", "TNF-alpha Signaling via NF-kB", "STAT1 (a)",
             "TP53 (a)", "Oxidative Phosphorylation", "DNA Repair", "Myc Targets V1")
    keepn = c("Ifn Gamma\nResponse", "TNFa\nvia NF-kB", "STAT1 (a)",
              "TP53 (a)", "Oxidative\nPhosphorylation", "DNA Repair", "Myc\nTargets V1")
    gdf = tsets %>%
        tidyr::gather("Gene set", "GSVA", -(sample:Aneuploidy)) %>%
        filter(!is.na(Type),
               `Gene set` %in% keep) %>%
        mutate(`Gene set` = factor(`Gene set`, levels=keep, ordered=TRUE))
    levels(gdf$`Gene set`) = keepn
    sdf = gdf %>% filter(!is.na(Driver))

    stars = c("<0.15"="⋅", "<0.05"="▪", "<0.001"="٭", "n.s."="×")
    a_v_b = function(x, y, a, b) {
        x1 = rep(NA, length(y))
        x1[x == a] = 0
        x1[x == b] = 1
        lm(y ~ x1) %>% broom::tidy() %>% filter(term == "x1") %>% pull(p.value)
    }
    get_p = function(df, var, t1, t2, t3, x="Gene set", y="GSVA") df %>%
        group_by(!!rlang::sym(x)) %>%
        summarize(l_vs_m = a_v_b(!!rlang::sym(var), !!rlang::sym(y), t1, t2),
                  m_vs_r = a_v_b(!!rlang::sym(var), !!rlang::sym(y), t2, t3)) %>%
        tidyr::gather("cmp", "p", -(!!rlang::sym(x))) %>%
        mutate(cmp = as.integer(factor(cmp, levels=c("l_vs_m", "m_vs_r"), ordered=TRUE)) - 1.5,
               !!rlang::sym(x) := as.integer(!!rlang::sym(x)) + sign(cmp)*0.15,
               sig = case_when(p<1e-3 ~ "<0.001", p<0.05 ~ "<0.05", p<0.15 ~ "<0.15", TRUE ~ "n.s."),
               sig = factor(sig, levels=names(stars)))
    sigs_type = get_p(gdf, "Type", "Myeloid", "B-like", "T-cell") %>% mutate(y = 0.7)
    aneup_type = get_p(tsets %>% mutate(x="1"), "Type", "Myeloid", "B-like", "T-cell", x="x", y="Aneuploidy")
    sigs_sub = get_p(sdf, "Driver", "Ebf1", "Ets1", "Erg") %>%
        mutate(y = ifelse(`Gene set` > 6.5, -0.35, 0.7))
    aneup_sub = get_p(subsets %>% mutate(x="1"), "Driver", "Ebf1", "Ets1", "Erg", x="x", y="Aneuploidy")

    common = function(map_bee) list(
        geom_boxplot(outlier.shape=NA),
        ggbeeswarm::geom_quasirandom(map_bee, color="black", alpha=0.3, dodge.width=.75, size=2.5, width=0.05),
        scale_shape_manual(values=c("Stat1"=25, "Pias1"=24), na.translate=TRUE, na.value=21, guide=guide_legend(order=2)),
        scale_discrete_manual("label", values=stars, drop=FALSE, name="p-value", guide=guide_legend(order=2)),
        theme(axis.title.x = element_blank())
    )
    set_common = function(sig_df, map_bee) c(
        list(geom_hline(yintercept=0, linetype="dashed", size=1, color="grey")),
        common(map_bee),
        list(guides(shape="none"),
             geom_text(data=sig_df, aes(x=`Gene set`, label=sig, y=y), size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE),
             plot_layout(tag_level="new"))
    )
    subvals = c("Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen1")
    p11 = ggplot(tsets, aes(x="Aneuploidy", y=Aneuploidy, fill=Type)) +
        common(aes(group=Type, shape=CIS)) +
        geom_text(data=aneup_type, aes(x=x, label=sig), y=0.57, size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
        scale_fill_discrete(guide=guide_legend(order=1)) +
        guides(shape="none")
    p21 = ggplot(subsets, aes(x="Aneuploidy", y=Aneuploidy, fill=Driver)) +
        common(aes(shape=CIS, group=Driver)) +
        geom_text(data=aneup_sub, aes(x=x, label=sig), y=0.4, size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
        scale_fill_manual(values=subvals, guide=guide_legend(order=1)) +
        guides(label="none")
    p12 = ggplot(gdf, aes(x=`Gene set`, y=GSVA, fill=Type)) +
        set_common(sigs_type, aes(shape=CIS, group=Type)) +
        scale_fill_discrete(guide=guide_legend(order=1))
    p22 = ggplot(sdf, aes(x=`Gene set`, y=GSVA, fill=Driver)) +
        set_common(sigs_sub, aes(shape=CIS, group=Driver)) +
        scale_fill_manual(values=subvals, guide=guide_legend(order=1)) +
        guides(label="none")
    top = p11 + p12 + plot_layout(widths=c(1,6.5), guides="collect")
    btm = p21 + p22 + plot_layout(widths=c(1,6.5), guides="collect")
    (top / btm)
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('c', 'cis', 'rds', '../cis_analysis/poisson.rds'),
        opt('e', 'expr', 'rds', '../expr_diff/de_Mad2PB.rds'),
        opt('g', 'gene_copies', 'rds', '../ploidy_compare/gene_copies.rds'),
        opt('h', 'gsva_hm', 'rds', '../data/gsva/mad2pb/MSigDB_Hallmark_2020.rds'),
        opt('d', 'gsva_dorothea', 'rds', '../data/gsva/mad2pb/DoRothEA.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig3-RNAseq.pdf')
    )

    aset = readRDS(args$aset)
    cis = readRDS(args$cis)$samples %>%
        filter(external_gene_name %in% c("Stat1", "Pias1")) %>%
        select(sample, CIS=external_gene_name)
    meta = aset$meta %>%
        left_join(cis) %>%
        transmute(sample=sample, Type=type, Subtype=subtype, CIS=CIS,
                  Aneuploidy=aneuploidy, Aneuploidy0.2=pmin(aneuploidy, 0.2))

    set_keep = c("Myc Targets V1", "Myc Targets V2", "Interferon Gamma Response",
                "Interferon Alpha Response", "Inflammatory Response", "DNA Repair",
                "Oxidative Phosphorylation", "TNF-alpha Signaling via NF-kB",
                "Mitotic Spindle", "TGF-beta Signaling", "STAT1 (a)", "TP53 (a)")
    gsva = t(narray::stack(readRDS(args$gsva_hm), readRDS(args$gsva_dorothea), along=1)) %>%
        as.data.frame() %>% tibble::rownames_to_column("sample") %>% as_tibble() %>%
        `[`(, c("sample", set_keep))

    gcs = readRDS(args$gene_copies) %>% t() %>%
        as.data.frame() %>% tibble::rownames_to_column("sample") %>% as_tibble() %>%
        select(sample, Myc_copies=`ENSMUSG00000022346`)
    gex = readRDS("../expr_diff/eset_Mad2PB.rds")$vs %>% t() %>%
        as.data.frame() %>% tibble::rownames_to_column("sample") %>% as_tibble() %>%
        select(sample, Myc_expr=Myc, Stat1_expr=Stat1, Pias1_expr=Pias1)

    dset = meta %>%
        dplyr::rename(Driver=Subtype) %>%
        left_join(gcs) %>%
        left_join(gex) %>%
        left_join(gsva)

    markers = readRDS("../expr_markers/markers.rds")
    diff_expr = readRDS(args$expr)

    #TODO:
    # quant cor black line +/- tissue adjustment
    # cor plot annotate insertions
    # volcano condition on STAT+Ifn (?)
    # volcano different colors for HMs, dorothea
    # myc copies -> myc targets? (maybe: does Myc targets assoc drop when conditioning on copies) [could do xy instead of @volc]

    umap = umap_types(markers)
    volc2 = aneup_volcano(diff_expr)
    gsa = set_tissue(dset) / gset_aneup(dset) + plot_layout(heights=c(1,1,2))

    asm = umap / (volc2 + gsa + plot_layout(widths=c(2,3))) +
        plot_layout(heights=c(1,2)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf(args$plotfile, 16.5, 15)
    print(asm)
    dev.off()
})
