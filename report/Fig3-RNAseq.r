library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
plt = import('plot')

umap_types = function(em) {
    set.seed(15208)
    meta = em$meta %>% filter(analysis_set & !is.na(type))
    vst = as.matrix(SummarizedExperiment::assay(em$vst[,meta$sample]))
    umap2 = uwot::umap(t(vst), n_components=2)
    dset = meta %>%
        mutate(umap1 = umap2[,1],
               umap2 = umap2[,2])
    dset$umap1[dset$umap1 > 0] = dset$umap1[dset$umap1 > 0] - 1.5 # empty space is meaningless, better vis
    dset$umap2[dset$umap2 > 0] = dset$umap2[dset$umap2 > 0] - 1.5

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
        ylab("adjusted p-value (FDR) Wald change") +
        theme(#axis.line.y = element_blank(),
              panel.grid.major = element_line(color="#efefef", size=0.5))
}

gset_aneup = function(gdf) {
    cols = c("Myeloid"="#f8766d", "T-cell"="#619cff", "Other"="#00ba38",
             "Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen3")
    gdf = gdf %>%
        mutate(Subtype = ifelse(is.na(Subtype), as.character(Type), as.character(Subtype))) %>%
        filter(Subtype != "B-like") # rare B-like
    p1 = ggplot(gdf, aes(x=`Interferon Gamma Response`, y=Aneuploidy0.2)) +
        geom_point(aes(color=Subtype, size=`STAT1 (a)`), alpha=0.5) +
        geom_smooth(aes(color=Subtype), method="lm", se=FALSE) +
        geom_smooth(method="lm", se=FALSE, color="black") +
        scale_color_manual(values=cols) +
        ylab("Aneuploidy") +
        plot_layout(tag_level="new")
    d1 = ggplot(gdf, aes(x=Type, y=`Interferon Gamma Response`)) +
        geom_boxplot(aes(fill=Type)) +
        coord_flip() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none")
    p2 = ggplot(gdf, aes(x=`Myc Targets V1`, y=Aneuploidy0.2)) +
        geom_point(aes(color=Subtype, size=`STAT1 (a)`), alpha=0.5) +
        geom_smooth(aes(color=Subtype), method="lm", se=FALSE) +
        geom_smooth(method="lm", se=FALSE, color="black") +
        theme(axis.title.y = element_blank()) +
        scale_color_manual(values=cols) +
        ylab("Aneuploidy") +
        plot_layout(tag_level="new")
    d2 = ggplot(gdf, aes(x=Type, y=`Myc Targets V1`)) +
        geom_boxplot(aes(fill=Type)) +
        coord_flip() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none")
    d3 = ggplot(gdf, aes(x=Type, y=Aneuploidy0.2)) +
        geom_boxplot(aes(fill=Type)) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none") +
        ylab("Aneuploidy") +
        plot_layout(tag_level="new")

    d1 + d2 + plot_spacer() + p1 + p2 + d3 +
        plot_layout(widths=c(5,5,1), heights=c(1,5), guides="collect")
}

set_tissue = function(sets) {
    tsets = sets %>% filter(!is.na(Type))
    subsets = sets %>% filter(!is.na(Subtype))
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
    sdf = gdf %>% filter(!is.na(Subtype))

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
    sigs_type = get_p(gdf, "Type", "Myeloid", "B-like", "T-cell")
    aneup_type = get_p(tsets %>% mutate(x="1"), "Type", "Myeloid", "B-like", "T-cell", x="x", y="Aneuploidy")
    sigs_sub = get_p(sdf, "Subtype", "Ebf1", "Ets1", "Erg") %>%
        mutate(y = ifelse(`Gene set` > 6.5, -0.35, 0.7))
    aneup_sub = get_p(subsets %>% mutate(x="1"), "Subtype", "Ebf1", "Ets1", "Erg", x="x", y="Aneuploidy")

    common = list(
        geom_boxplot(outlier.shape=NA),
        ggbeeswarm::geom_quasirandom(color="black", alpha=0.3, shape=21, dodge.width=.75, size=2.5, width=0.05),
        scale_fill_discrete(guide=guide_legend(order=1)),
        scale_discrete_manual("label", values=stars, drop=FALSE, name="p-value", guide=guide_legend(order=2)),
        theme(axis.title.x = element_blank())
    )
    p1 = ggplot(tsets, aes(x="Aneuploidy", y=Aneuploidy, fill=Type)) +
        common +
        geom_text(data=aneup_type, aes(x=x, label=sig), y=0.57, size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE)
    p2 = ggplot(gdf, aes(x=`Gene set`, y=GSVA, fill=Type)) +
        geom_hline(yintercept=0, linetype="dashed", size=1, color="grey") +
        common +
        geom_text(data=sigs_type, aes(x=`Gene set`, label=sig), y=0.7, size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
        plot_layout(tag_level="new")
    p3 = ggplot(subsets, aes(x="Aneuploidy", y=Aneuploidy, fill=Subtype)) +
        common +
        geom_text(data=aneup_sub, aes(x=x, label=sig), y=0.4, size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
        scale_fill_manual(values=c("Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen3")) +
        guides(label=FALSE)
    p4 = ggplot(sdf, aes(x=`Gene set`, y=GSVA, fill=Subtype)) +
        geom_hline(yintercept=0, linetype="dashed", size=1, color="grey") +
        common +
        geom_text(data=sigs_sub, aes(x=`Gene set`, label=sig, y=y), size=5, hjust=0.5, vjust=0.5, inherit.aes=FALSE) +
        scale_fill_manual(values=c("Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen3")) +
        guides(label=FALSE) +
        plot_layout(tag_level="new")
    top = p1 + p2 + plot_layout(widths=c(1,6.5), guides="collect")
    btm = p3 + p4 + plot_layout(widths=c(1,6.5), guides="collect")
    (top / btm)
}

ifn_myc_condition = function(sets, ghm, gdo) {
    narray::intersect(sets$sample, ghm, gdo, along=2)
    gsva = t(rbind(ghm, gdo))

    #todo: (1) aneup ~ Stat1_Ifn + gsva_cats
    #      (2) aneup ~ Myc_copies + gsva_cats
    # (3) plot wald stats on x/y, show that Ifn does not explain myc but copies do more
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
    meta = aset$meta

    gcs = readRDS(args$gene_copies)
    gex = readRDS("../expr_diff/eset_Mad2PB.rds")$vs

    ghm = readRDS(args$gsva_hm)
    gdo = readRDS(args$gsva_dorothea)
    narray::intersect(ghm, gdo, along=2)
    sets = cbind.data.frame(
        sample = colnames(ghm),
        Type = meta$type[match(colnames(ghm), meta$sample)],
        Subtype = meta$subtype[match(colnames(ghm), meta$sample)],
        Myc_expr = gex["Myc", match(colnames(ghm), colnames(gex))],
        Myc_copies = pmin(gcs["ENSMUSG00000022346", match(colnames(ghm), colnames(gcs))], 3),
        Aneuploidy = meta$aneuploidy[match(colnames(ghm), meta$sample)],
        Aneuploidy0.2 = pmin(meta$aneuploidy[match(colnames(ghm), meta$sample)], 0.2),
        t(ghm[c("Myc Targets V1", "Myc Targets V2", "Interferon Gamma Response",
                "Interferon Alpha Response", "Inflammatory Response", "DNA Repair",
                "Oxidative Phosphorylation", "TNF-alpha Signaling via NF-kB",
                "Mitotic Spindle", "TGF-beta Signaling"),]),
        t(gdo[c("STAT1 (a)", "TP53 (a)"),])
    ) %>% as_tibble()

    #todo: move the subtype annotations to the actual metadata
    cis = readRDS(args$cis)$samples # todo: add stat1 ins to lm plot

    markers = readRDS("../expr_markers/markers.rds")
    diff_expr = readRDS(args$expr)

    #TODO:
    # quant cor black line +/- tissue adjustment
    # cor plot annotate insertions
    # volcano condition on STAT+Ifn (?)
    # volcano different colors for HMs, dorothea
    # switch Ighm for monocyte marker?
    # myc copies -> myc targets? (maybe: does Myc targets assoc drop when conditioning on copies) [could do xy instead of @volc]

    umap = umap_types(markers)
    volc2 = aneup_volcano(diff_expr)
    gsa = set_tissue(sets) / gset_aneup(sets) + plot_layout(heights=c(1,1,2))

    asm = umap / (volc2 + gsa + plot_layout(widths=c(2,3))) +
        plot_layout(heights=c(1,2)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf(args$plotfile, 16.5, 15)
    print(asm)
    dev.off()
})
