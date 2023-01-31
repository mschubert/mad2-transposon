library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')
gset = import('genesets')
tcga = import('data/tcga')
sys = import('sys')

get_muts = function(aneup, ex) {
    load_fl = function(coh) tcga$mutations(coh) %>%
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol, vclass=Variant_Classification)
    m = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows() %>%
        inner_join(aneup) %>%
        group_by(aneup_class, Sample) %>%
            mutate(tot = n_distinct(gene)) %>%
        ungroup() %>%
        filter(tot <= 5000) %>%
        inner_join(ex)
}

test_fet = function(muts) {
    genes = table(muts$gene)
    genes = names(genes)[genes > 50]

    test_one = function(g) {
        df = muts %>% filter(g == gene) %>% group_by(aneup_class, Sample, glen) %>%
            summarize(.gene = -log10(tot[1]))
        broom::tidy(lm(.gene ~ aneup_class, data=df)) %>% mutate(glen=df$glen[1])
    }
    sapply(genes, test_one, simplify=FALSE) %>% bind_rows(.id="gene") %>%
        filter(term == "aneup_class(0.1,Inf]") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>%
        filter(std.error < 0.2)
}

distr_plot = function(res) {
#    res$circle = res$gene %in% c("JAK1", "STAT1", "B2M", "TP53", "TTN", "CDKN2A", "HLA-A", "HLA-B", "ERBB2", "KRAS", "MYC")
    res$circle = res$gene %in% c("JAK1", "STAT1", "TP53", "TTN")
#    plt$volcano(res, ceil=1e-15, label_top=50)

    #TODO: figure out why TP53 vs. CDKN2A so different here
    res2 = res[res$circle,]
    ggplot(res, aes(x=estimate)) +
        geom_density() +
        geom_point(data=res2, y=0) +
        ggrepel::geom_text_repel(data=res2, y=0, aes(label=gene)) #+
#        coord_cartesian(xlim=c(-2e-4, 7e-4))
}

volcano_hallmarks = function(res) {
    sets_hm = gset$get_human("MSigDB_Hallmark_2020")
    hits_hm = res$gene[res$estimate>quantile(res$estimate,0.7)]

    fet_hm = gset$test_fet(hits_hm, sets_hm, valid=res$gene)
    lm_hm = gset$test_lm(res, stat="estimate", sets_hm)

    plt$volcano(fet_hm %>% mutate(estimate=log2(estimate))) + ggtitle("FET") |
    plt$volcano(lm_hm) + ggtitle("lm")
}

volcano_go = function(res) {
    sets_go = gset$get_human("GO_Biological_Process_Tree")
    hits_go = res$gene[res$estimate>quantile(res$estimate,0.7)]

    fet_go = gset$test_fet(hits_go, sets_go, valid=res$gene, min=4)
    lm_go = gset$test_lm(res, stat="estimate", gset$filter(sets_go, max=250))

    plt$volcano(fet_go %>% mutate(estimate=log2(estimate))) + ggtitle("FET") |
    plt$volcano(lm_go) + ggtitle("lm")
}

sys$run({
    aneup = readRDS("mut_frac.rds")$aneup
    ex = readRDS("canonical_exon_length.rds")

    muts = get_muts(aneup, ex)
    res = test_fet(muts)

    saveRDS(list(muts=muts, res=res), file="mut_distr.rds")
})
