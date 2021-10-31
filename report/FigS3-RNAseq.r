library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
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

set_cor = function(sets) {
    smat = data.matrix(sets[-(1:2)])
    rownames(smat) = sets$sample
    meta2 = meta[match(sets$sample, meta$sample),]
    smat = cbind(smat, narray::mask(meta2$type))

    cc = cor(smat)
    ggcorrplot::ggcorrplot(cc, hc.order=TRUE, method="circle")
}

EtsErg_subtype = function() {
    dset = readRDS("../expr_diff/eset_Mad2PB.rds")
    meta = as.data.frame(SummarizedExperiment::colData(dset$eset))
    idx = cbind(meta, t(dset$vs[c("Ets1", "Erg"),]))
    idx = idx[idx$type != "unknown",]

    p1 = ggplot(idx, aes(x=Erg, y=Ets1)) +
        geom_point(aes(size=aneuploidy, fill=type), shape=21) +
        ggrepel::geom_text_repel(aes(label=sample), size=2) +
        labs(title = "Mouse cohort expression",
             fill = "Mouse tumor type",
             size = "Aneuploidy")

    hutype = c(
        "T-ALL" = "#2b8cbe",
        "ALL with hyperdiploid karyotype" = "#54278f",
        "c-ALL/Pre-B-ALL with t(9;22)" = "#d9f0a3",
        "c-ALL/Pre-B-ALL without t(9;22)" = "#78c679",
        "Pro-B-ALL with t(11q23)/MLL" = "#238443",
        "mature B-ALL with t(8;14)" = "#ffffe5",
        "AML with normal karyotype + other abnormalities" = "#7f2704",
        "AML with complex aberrant karyotype" = "#feb24c"
    )

    dset2 = readRDS("../data/arrayexpress/E-GEOD-13159.rds")
    expr = Biobase::exprs(dset2)
    rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
    meta = Biobase::pData(dset2) %>%
        transmute(sample = Array.Data.File,
                  type = FactorValue..LEUKEMIA.CLASS.)
    idx2 = cbind(meta, t(expr[c("ETS1", "ERG"),]))
    idx2 = idx2[idx2$type %in% names(hutype),]

    p2 = ggplot(idx2, aes(x=ERG, y=ETS1)) +
        geom_point(aes(fill=type, size=aneuploidy), alpha=0.5, shape=21, size=3) + # size=aneuploidy
        labs(title = "Human leukemia cohort (MILE)",
             fill = "Human leukemia type",
             size = "Aneuploidy") +
        scale_fill_manual(values=hutype)

    p1 + p2 + plot_layout(guides="collect")
}

#fixme: defunct
EtsErg_TPS = function() {
    genes = c("Ets1", "Erg", "Stat1", "Pias1")
    dset = readRDS("../expr_diff/eset_Mad2PB.rds")

    meta = as.data.frame(SummarizedExperiment::colData(dset$eset))

    expr = reshape2::melt(dset$vs[genes,]) %>%
        dplyr::rename(gene=Var1, sample=Var2, expr=value) %>%
        inner_join(meta %>% select(sample, type)) %>%
        filter(type != "unknown") %>%
        mutate(type = as.character(type))

    eplot = ggplot(expr, aes(x=type, y=expr)) +
        geom_boxplot(aes(fill=type), outlier.shape=NA) +
        facet_wrap(~ gene, scales="free_y", ncol=1)

    #FIXME: annoying empty plots
    c1 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "Myeloid"])), main="Myeloid")
    c2 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "Other"])), main="B-like")
    c3 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "T-cell"])), main="T-ALL")
    { c1 + c2 + c3 + plot_layout(ncol=1) } + { eplot } + plot_layout(nrow=1)
}

#fixme: defunct
EtsErg_MILE = function() {
    mad2pb = io$load("../expr_diff/eset_Mad2PB.RData")
    annot = as.data.frame(SummarizedExperiment::colData(mad2pb$eset))

    mile = io$load("../expr_diff/eset_MILE.RData")#$expr # is there a coarse def?
    meta = mile$meta

    genes = c("Ets1", "Erg", "Stat1", "Pias1", "Ifng")#, "Stat3", "Pten", "Notch1")
    c1 = cor(t(mad2pb$vs[genes, annot$type=="Myeloid"]))
    c2 = cor(t(mad2pb$vs[genes, annot$type=="B-like"]))
    c3 = cor(t(mad2pb$vs[genes, annot$type=="T-cell"]))

    genes = c("ETS1", "ERG", "STAT1", "PIAS1", "Ifng")#, "STAT3", "PTEN", "NOTCH1")
    m1 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="Myeloid"]))
    m2 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="B_like"]))
    m3 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="T_ALL"]))

    corrplot::corrplot(c1, tl.cex=2, tl.col="black")
    corrplot::corrplot(c2, tl.cex=2, tl.col="black")
    corrplot::corrplot(c3, tl.cex=2, tl.col="black")

    corrplot::corrplot(m1, tl.cex=2, tl.col="black")
    corrplot::corrplot(m2, tl.cex=2, tl.col="black")
    corrplot::corrplot(m3, tl.cex=2, tl.col="black")
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('c', 'cis', 'rds', '../cis_analysis/poisson.rds'),
        opt('e', 'expr', 'rds', '../expr_diff/de_Mad2PB.rds'),
        opt('g', 'gene_copies', 'rds', '../ploidy_compare/gene_copies.rds'),
        opt('h', 'gsva_hm', 'rds', '../data/gsva/mad2pb/MSigDB_Hallmark_2020.rds'),
        opt('d', 'gsva_dorothea', 'rds', '../data/gsva/mad2pb/DoRothEA.rds'),
        opt('p', 'plotfile', 'pdf', 'FigS3-RNAseq.pdf')
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
        type = meta$type[match(colnames(ghm), meta$sample)],
        subtype = meta$subtype[match(colnames(ghm), meta$sample)],
        Myc_expr = gex["Myc", match(colnames(ghm), colnames(gex))],
        Myc_copies = pmin(gcs["ENSMUSG00000022346", match(colnames(ghm), colnames(gcs))], 3),
        Aneuploidy = pmin(meta$aneuploidy[match(colnames(ghm), meta$sample)], 0.2),
        t(ghm[c("Myc Targets V1", "Myc Targets V2", "Interferon Gamma Response",
                "Interferon Alpha Response", "Inflammatory Response", "DNA Repair",
                "Oxidative Phosphorylation", "TNF-alpha Signaling via NF-kB",
                "Mitotic Spindle", "TGF-beta Signaling"),]),
        t(gdo[c("STAT1 (a)", "TP53 (a)"),])
    ) %>% as_tibble() %>% na.omit()

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

    r1 = wrap_elements(marker_pca(markers))
    r2 = wrap_elements(EtsErg_subtype())

    asm = (r1 / r2) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 12)
    print(asm)
    dev.off()
})
