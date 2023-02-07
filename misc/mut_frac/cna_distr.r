library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')
gset = import('genesets')
tcga = import('data/tcga')
sys = import('sys')

get_cnas = function(aneup) {
    cnas = tcga$cna_gistic(thresh=TRUE) %>%
        reshape2::melt() %>%
        dplyr::rename(Sample=Var2, gene=Var1, gistic=value) %>%
        as_tibble() %>%
        inner_join(aneup) %>%
        filter(substr(Sample, 14, 16) == "01A") %>%
        group_by(aneup_class, Sample) %>%
            mutate(tot = n_distinct(gene[abs(gistic) == 2])) %>%
        ungroup()
}

test_fet = function(cnas) {
    cnas2 = cnas %>%
        filter(gistic == -2) %>%
        group_by(gene) %>%
            filter(n_distinct(aneup_class) == 2,
                   n_distinct(Sample) >= 50) %>%
        ungroup() %>%
        mutate(gene = droplevels(gene)) %>%
        split(.$gene)

    test_one = function(df) {
        df2 = df %>% mutate(.gene = -log10(tot))
        broom::tidy(lm(.gene ~ aneup_class, data=df2))
    }
    lapply(cnas2, test_one) %>% bind_rows(.id="gene") %>%
        filter(term == "aneup_class(0.1,Inf]") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>%
        filter(std.error < 0.2)
}

distr_plot = function(res) {
    res$circle = res$gene %in% c("JAK1", "STAT1", "TP53", "TTN", "CDKN2A", "IFNA1", "IFNB1")
#    plt$volcano(res, ceil=1e-15, label_top=50)

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
    plot = plt$volcano(fet_hm %>% mutate(estimate=log2(estimate))) + ggtitle("FET") |
        plt$volcano(lm_hm) + ggtitle("lm")

    list(fet=fet_hm, lm=lm_hm, plot=plot)
}

volcano_go = function(res) {
    sets_go = gset$get_human("GO_Biological_Process_Tree")
    hits_go = res$gene[res$estimate>quantile(res$estimate,0.7)]

    fet_go = gset$test_fet(hits_go, sets_go, valid=res$gene, min=4)
    lm_go = gset$test_lm(res, stat="estimate", gset$filter(sets_go, max=250))
    plot = plt$volcano(fet_go %>% mutate(estimate=log2(estimate))) + ggtitle("FET") |
        plt$volcano(lm_go) + ggtitle("lm")

    list(fet=fet_go, lm=lm_go, plot=plot)
}

sys$run({
    aneup = readRDS("mut_frac.rds")$aneup

    cnas = get_cnas(aneup)
    res = test_fet(cnas)

    vhm = volcano_hallmarks(res)
    vgo = volcano_go(res)

    pdf("cna_distr.pdf", 14, 8)
    print(distr_plot(res))
    print(vhm$plot)
    print(vgo$plot)
    dev.off()

    saveRDS(list(cnas=cnas, res=res, fet_hm=vhm$fet), file="cna_distr.rds")
})
