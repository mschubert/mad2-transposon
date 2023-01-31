library(dplyr)
library(ggplot2)
library(patchwork)
tcga = import('data/tcga')
sys = import('sys')

all_muts_volc = function(aneup, ex) {
    plt = import('plot')
    gset = import('genesets')

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
    genes = table(m$gene)
    genes = names(genes)[genes > 50]

    test_one = function(g) {
        df = m %>% filter(g == gene) %>% group_by(aneup_class, Sample, glen) %>% summarize(.gene = -log10(tot[1]))
        broom::tidy(lm(.gene ~ aneup_class, data=df)) %>% mutate(glen=df$glen[1])
    }
    res = sapply(genes, test_one, simplify=FALSE) %>% bind_rows(.id="gene") %>%
        filter(term == "aneup_class(0.1,Inf]") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>%
        filter(std.error < 0.2)
#    res$circle = res$gene %in% c("JAK1", "STAT1", "B2M", "TP53", "TTN", "CDKN2A", "HLA-A", "HLA-B", "ERBB2", "KRAS", "MYC")
    res$circle = res$gene %in% c("JAK1", "STAT1", "TP53", "TTN")
    plt$volcano(res, ceil=1e-15, label_top=50)

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "GO_Biological_Process_Tree"))
    f1 = gset$test_fet(res$gene[res$estimate>quantile(res$estimate,0.7)], sets[[1]], valid=res$gene)
    plt$volcano(f1 %>% mutate(estimate=log2(estimate)))
    s1 = gset$test_lm(res, stat="estimate", sets[[1]])
    plt$volcano(s1)
    f2 = gset$test_fet(res$gene[res$estimate>quantile(res$estimate,0.7)], sets[[2]], valid=res$gene, min=4)
    plt$volcano(f2 %>% mutate(estimate=log2(estimate)))
    s2 = gset$test_lm(res, stat="est2", gset$filter(sets[[2]], max=250))
    plt$volcano(s2)

    #TODO: figure out why TP53 vs. CDKN2A so different here
    res2 = res[res$circle,]
    ggplot(res, aes(x=estimate)) +
        geom_density() +
        geom_point(data=res2, y=0) +
        ggrepel::geom_text_repel(data=res2, y=0, aes(label=gene)) #+
#        coord_cartesian(xlim=c(-2e-4, 7e-4))
}

sys$run({
    aneup = readRDS("mut_frac.rds")$aneup
    ex = readRDS("canonical_exon_length.rds")

    all_muts_volc(aneup, ex)
})
