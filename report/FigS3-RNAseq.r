library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')

marker_pca = function(markers) {
    one_pca = function(mm, smps) {
        plt$pca(prcomp(mm[smps,]), aes(x=PC1, y=PC2), annot=meta[smps,],
                biplot=TRUE, bi_size=2.5, bi_color="darkviolet", bi_arrow=0.2) +
            geom_point(aes(color=type, size=aneuploidy), alpha=0.8) +
            ggrepel::geom_text_repel(aes(label=sample), size=2, segment.alpha=0.3) +
            scale_color_discrete(drop=FALSE) +
            theme_classic()
    }
    meta = markers$meta %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.3))
    mmat = t(SummarizedExperiment::assay(markers$vst))

    # separate T-cells
    mm1 = mmat[,c("Cd3g", "Cd3e", "Rag1", "Sox4", "Dntt", "Xrcc6", "Notch1", "Ly6a",
                  "Myc", "Cd4", "Cd8a", "Il7", "Cd34")]
    p1 = one_pca(mm1, TRUE)

    # separate myeloid from B-like
    smps = is.na(meta$type) | meta$type != "T-cell"
    mm2 = mmat[,c("Ebf1", "Myc", "Ly6a", "Ly6g", "Cd14", "Nkg7", "Ms4a1", "Cd34", "Lyz1", "Ms4a7")]
    p2 = one_pca(mm2, smps)

    # separate Ets-Erg-Myc-Cd19 axis
    smps = is.na(meta$type) | meta$type %in% "B-like"
    mm3 = mmat[,c("Ebf1", "Ets1", "Erg", "Myc", "Stat1", "Kit", "Tmem184a", "Cd34", "Pax5", "Cd19")]
    p3 = one_pca(mm3, smps)

    p1 + (p2 + plot_layout(tag_level="new")) +
        (p3 + plot_layout(tag_level="new")) +
        plot_layout(nrow=1, guides="collect")
}

EtsErg_subtype = function(tps, mile, hutype) {
    p1 = ggplot(tps, aes(x=Erg, y=Ets1)) +
        geom_point(aes(size=aneuploidy, fill=type), shape=21) +
        ggrepel::geom_text_repel(aes(label=sample), size=2) +
        labs(title = "Mouse cohort expression",
             fill = "Mouse tumor type",
             size = "Aneuploidy")

    p2 = ggplot(mile, aes(x=ERG, y=ETS1)) +
        geom_point(aes(fill=type, size=aneuploidy), alpha=0.5, shape=21) +
        labs(title = "Human leukemia cohort (MILE)",
             fill = "Human leukemia type",
             size = "Aneuploidy") +
        scale_fill_manual(values=hutype)

    p1 + p2 + plot_layout(guides="collect")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'markers', 'rds', '../expr_markers/markers.rds'),
        opt('e', 'eset', 'rds', '../expr_diff/eset_Mad2PB.rds'),
        opt('h', 'human', 'rds', '../data/arrayexpress/E-GEOD-13159.rds'),
        opt('p', 'plotfile', 'pdf', 'FigS3-RNAseq.pdf')
    )

    cols = c("Myeloid"="#f8766d", "T-cell"="#619cff", "B-like"="#00ba38",
             "Ets1"="chartreuse3", "Erg"="forestgreen", "Ebf1"="darkolivegreen3")
    markers = readRDS(args$markers)
    aset = readRDS("../ploidy_compare/analysis_set.rds")$meta
    dset = readRDS(args$eset)
    meta = as.data.frame(SummarizedExperiment::colData(dset$eset))
    tps = cbind(meta, t(dset$vs[c("Ets1", "Erg"),]))
    tps = tps[tps$type != "unknown",] %>%
        inner_join(aset %>% select(sample, subtype)) %>%
        mutate(subtype = ifelse(type=="Other", as.character(subtype), as.character(type))) %>%
        filter(!is.na(type))

    hutype = c(
        "T-ALL" = "#2b8cbe",
        "ALL with hyperdiploid karyotype" = "#54278f",
        "c-ALL/Pre-B-ALL with t(9;22)" = "#d9f0a3",
        "c-ALL/Pre-B-ALL without t(9;22)" = "#78c679",
        "Pro-B-ALL with t(11q23)/MLL" = "#238443",
        "mature B-ALL with t(8;14)" = "#ffffe5",
        "AML with normal karyotype" = "#7f2704",
        "AML with complex aberrant karyotype" = "#feb24c"
    )
    dset2 = readRDS(args$human)
    expr = Biobase::exprs(dset2)
    rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
    meta = Biobase::pData(dset2) %>%
        transmute(sample = Array.Data.File,
                  type = FactorValue..LEUKEMIA.CLASS.)
    mile = cbind(meta, t(expr[c("ETS1", "ERG"),]))
    mile$type = sub(" + other abnormalities", "", mile$type, fixed=TRUE)
    mile = mile[mile$type %in% names(hutype),] %>%
        inner_join(import('io')$load("../misc/subtype_overlay/aneup_scores_mile.RData")$aneup)
#        inner_join(readRDS("../misc/subtype_overlay/aneup_scores_mile.rds")$aneup)

    aneups = list(`Mouse transposon cohort` = tps %>% select(sample, type=subtype, aneuploidy),
                  `Human leukemia cohort (MILE)` = mile %>% select(sample, type, aneuploidy)) %>%
        bind_rows(.id="dset") %>%
        filter(!is.na(type)) %>%
        as_tibble()
    pa = ggplot(aneups, aes(x=forcats::fct_reorder(type, aneuploidy), y=aneuploidy, fill=type)) +
        geom_boxplot() +
        facet_grid(. ~ dset, scales="free", space="free") +
        scale_x_discrete(guide = guide_axis(n.dodge=2)) +
        scale_fill_manual(values=c(cols, hutype), guide="none") +
        theme(axis.title.x = element_blank()) +
        ylab("Aneuploidy")

    r1 = wrap_elements(marker_pca(markers))
    r2 = wrap_elements(EtsErg_subtype(tps, mile, hutype))

    asm = (r1 / r2 / wrap_elements(pa)) +
        plot_layout(heights=c(1,1.2,0.8)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 16)
    print(asm)
    dev.off()
})
