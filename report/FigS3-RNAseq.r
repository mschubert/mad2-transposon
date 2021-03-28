library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
gset = import('genesets')
plt = import('plot')
gnet = import('tools/genenet')

marker_pca = function() {
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

    pdf(args$plotfile, 16, 15)
#    print(asm)
    dev.off()
})
