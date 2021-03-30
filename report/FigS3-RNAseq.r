library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
plt = import('plot')
gnet = import('tools/genenet')

marker_pca = function(markers) {
    one_pca = function(mm, smps) {
        plt$pca(prcomp(mm[smps,]), aes(x=PC1, y=PC2), annot=meta[smps,],
                biplot=TRUE, bi_size=3.5, bi_color="black") +
            geom_point(aes(color=type, size=aneuploidy), alpha=0.8) +
            ggrepel::geom_text_repel(aes(label=sample), size=2) +
            scale_color_discrete(drop=FALSE) +
            theme_classic()
    }
    meta = markers$meta %>%
        mutate(type = factor(type),
               aneuploidy = pmin(aneuploidy, 0.3))
    mmat = t(SummarizedExperiment::assay(markers$vst))

    # separate T-cells
    mm1 = mmat[,c("Cd3g", "Cd3e", "Rag1", "Sox4", "Dntt", "Xrcc6", "Notch1", "Ly6a",
                  "Myc", "Cd4", "Cd8a", "Il7", "Cd34")]
    p1 = one_pca(mm1, TRUE)

    # separate myeloid from B-like
    smps = !is.na(meta$type) & meta$type != "T-cell"
    mm2 = mmat[,c("Ebf1", "Myc", "Ly6a", "Ly6g", "Cd14", "Nkg7", "Ms4a1", "Cd34", "Lyz1", "Ms4a7")]
    p2 = one_pca(mm2, smps)

    # separate Ets-Erg-Myc-Cd19 axis
    smps = !is.na(meta$type) & meta$type %in% "Other"
    mm3 = mmat[,c("Ebf1", "Ets1", "Erg", "Myc", "Stat1", "Kit", "Tmem184a", "Cd34", "Pax5", "Cd19")]
    p3 = one_pca(mm3, smps)

    p1 + p2 + p3 + plot_layout(nrow=1, guides="collect")
}

set_cor = function(sets) {
    smat = data.matrix(sets[-(1:2)])
    rownames(smat) = sets$sample
    meta2 = meta[match(sets$sample, meta$sample),]
    smat = cbind(smat, narray::mask(meta2$type))

    cc = cor(smat)
    ggcorrplot::ggcorrplot(cc, hc.order=TRUE, method="circle")
}

cor_net = function(sets, cis) {
#TODO;
# should we do
# (1) selected HMs boxplot between tissue -> t-cells have a lot already (B: depends Erg; myeloid: not)
# (2) GEX correlation + ins changing expr of category

    types = narray::mask(sets$type) + 0
    smat = cbind(types, data.matrix(sets[-c(1,2)]))
    rownames(smat) = rownames(types) = sets$sample

    test_cis = c("Ets1", "Erg", "Trp53", "Stat1", "Pten", "Ikzf1", "Crebbp")
    cis2 = cis %>%
        filter(external_gene_name %in% test_cis) %>%
        mutate(has_ins = 1) %>%
        narray::construct(has_ins ~ sample + external_gene_name, fill=0)

    narray::intersect(sets$sample, types, smat, cis2, along=1)

    stype = sets$type
    res = st$lm(smat ~ stype + cis2) %>%
        as_tibble() %>%
        filter(term == "cis2") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    ggplot(res, aes(x=smat, y=cis2, fill=estimate)) +
        geom_tile() +
        scale_fill_distiller(palette="RdBu", lim=c(-0.5,0.5)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    cmat = cbind(types, cis2)
    smat2 = data.matrix(sets[-c(1,2)])
    res3 = st$lm(smat2 ~ 0 + cmat, atomic="cmat") %>%
        as_tibble() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    ggplot(res3, aes(x=smat2, y=term, fill=estimate)) +
        geom_tile() +
        scale_fill_distiller(palette="RdBu", lim=c(-1,1)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    res2 = st$lm(smat ~ types) %>%
        as_tibble() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    x = res2 %>% filter(! smat %in% c("Other", "T-cell", "Myeloid"))
    ggplot(x, aes(x=smat, y=types, fill=estimate)) +
        geom_tile() +
        scale_fill_distiller(palette="RdBu", lim=c(-1,1)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    #todo: add myc CN, expr?

#    smat = data.matrix(sets[-(1:2)])
#    smat = smat[,-c(3,5,6)]

    smat5 = smat[,-c(1,2,3,8,10,11,15,16,5)] # 16=stat1, 14=mit spind, 16 tgfb, 6 aneup, 5 myc copies

    gnet$plot_bootstrapped_pcor(smat5, fdr=0.5, n=500, show_edge_if=50, layout="stress")
}

#fixme: defunct
EtsErg_TPS = function() {
    genes = c("Ets1", "Erg", "Stat1", "Pias1")
    dset = io$load("../expr_diff/eset_Mad2PB.RData")

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
    #{ c1 + c2 + c3 + plot_layout(ncol=1) } + { eplot } + plot_layout(nrow=1)
}

#fixme: defunct
EtsErg_MILE = function() {
    mad2pb = io$load("../expr_diff/eset_Mad2PB.RData")
    annot = as.data.frame(SummarizedExperiment::colData(mad2pb$eset))

    mile = io$load("../expr_diff/eset_MILE.RData")#$expr # is there a coarse def?
    meta = mile$meta

    genes = c("Ets1", "Erg", "Stat1", "Pias1", "Ifng")#, "Stat3", "Pten", "Notch1")
    c1 = cor(t(mad2pb$vs[genes, annot$type=="Myeloid"]))
    c2 = cor(t(mad2pb$vs[genes, annot$type=="Other"]))
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
        Myc_expr = gex["Myc", match(colnames(ghm), colnames(gex))],
        Myc_copies = pmin(gcs["ENSMUSG00000022346", match(colnames(ghm), colnames(gcs))], 3),
        Aneuploidy = pmin(meta$aneuploidy[match(colnames(ghm), meta$sample)], 0.2),
        t(ghm[c("Myc Targets V1", "Myc Targets V2", "Interferon Gamma Response",
                "Interferon Alpha Response", "Inflammatory Response", "DNA Repair",
                "Oxidative Phosphorylation", "TNF-alpha Signaling via NF-kB",
                "Mitotic Spindle", "TGF-beta Signaling"),]),
        t(gdo[c("STAT1 (a)", "TP53 (a)"),])
    ) %>% as_tibble() %>% na.omit()

    #todo: move the subtype annotations to the actual metadata
    gdf = sets %>%
        mutate(subtype = case_when(
            sample %in% c("157s", "404s", "411s", "412s", "416s", "424s", "425s", "428s",
                          "432s", "437s", "461s", "476s", "482s", "627s", "632s") ~ "Ets1",
            sample %in% c("402s", "405s", "409s", "413s", "417s", "422s", "429s", "431s",
                          "434s", "435s", "442s") ~ "Erg",
            sample %in% c("421s", "467s", "485s", "609s", "613t", "620s", "622s") ~ "Ebf1",
            TRUE ~ NA_character_
        )) %>%
        mutate(subtype = factor(subtype, levels=c("Ebf1", "Ets1", "Erg"))) %>%
        select(sample, subtype, everything())

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

    asm = marker_pca(markers) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 5)
    print(asm)
    dev.off()
})
