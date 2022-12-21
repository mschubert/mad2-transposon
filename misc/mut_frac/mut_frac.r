library(dplyr)
library(ggplot2)
library(patchwork)
tcga = import('data/tcga')
sys = import('sys')

get_aneup_purity = function() {
    aneup = tcga$aneuploidy() %>%
        select(Sample, aneup=aneup_log2seg) %>%
        inner_join(tcga$purity() %>% select(Sample, cohort, purity=estimate)) %>%
        filter(substr(Sample, 14, 16) == "01A",
               !is.na(aneup) & !is.na(purity),
               purity > 0.75) %>%
        mutate(aneup = aneup / purity,
               aneup_class = cut(aneup, c(0,0.1,Inf)))
}

get_cnas = function(aneup) {
    cnas = tcga$cna_gistic(thresh=TRUE) %>%
        reshape2::melt() %>%
        dplyr::rename(Sample=Var2, gene=Var1, gistic=value) %>%
        as_tibble() %>%
        filter(substr(Sample, 14, 16) == "01A")
    cfrac = cnas %>%
        group_by(Sample) %>%
            summarize(n_eup = sum(gistic == 1),
                      n_amp = sum(gistic == 2),
                      n_del = sum(gistic == -2),
                      `TP53 loss` = as.integer(gistic[gene == "TP53"] == -2),
    #                  `MYC amp` = as.integer(gistic[gene == "MYC"] == 2),
                      `IFNA1 loss` = as.integer(gistic[gene == "IFNA1"] == -2)) %>%
    #                  `CDKN2A loss` = as.integer(gistic[gene == "CDKN2A"] == -2)) %>%
        ungroup()
    xx = inner_join(aneup, cfrac) %>%
        tidyr::gather("gene", "cna", -(Sample:n_del)) %>%
        mutate(gene = factor(gene, levels=c("MYC amp", "IFNA1 loss", "TP53 loss")))
}

plot_cnas = function(xx) {
    p1 = ggplot(xx, aes(x=factor(cna), y=aneup, color=factor(cna))) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_color_brewer(palette="Set1", name="Genomic\nevent") +
        labs(x = "Presence of event",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("0", "1")))

    p2 = ggplot(xx, aes(x=aneup_class, y=cna/(n_del), color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_y_log10() +
        labs(x = "Aneuploidy",
             y = "As fraction of genes lost") +
        scale_color_manual(values=setNames(c("#cc9933", "#5500aa"), c("(0,0.1]", "(0.1,Inf]")), name="Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("(0,0.1]", "(0.1,Inf]")))

    p1 | p2
}

get_muts = function(aneup) {
    load_fl = function(coh) tcga$mutations(coh) %>%
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol, vclass=Variant_Classification)
    m = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows() %>%
        group_by(cohort, Sample) %>%
            mutate(n = n()) %>%
            summarize(tot = unique(n),
                      #STAT1 = as.integer("STAT1" %in% gene[vclass != "Silent"]),
                      TP53 = as.integer("TP53" %in% gene[vclass != "Silent"]),
                      TTN = as.integer("TTN" %in% gene[vclass != "Silent"])) %>%
        ungroup()
    both = aneup %>%
        inner_join(m) %>%
        tidyr::gather("gene", "mut", -(Sample:tot)) %>%
        filter(tot <= 5000)
}

plot_muts = function(both) {
    p1 = ggplot(both, aes(x=factor(mut), y=aneup, color=factor(mut))) +
        geom_boxplot(outlier.shape=NA) +
        scale_color_brewer(palette="Set1", name="Genomic\nevent") +
        facet_wrap(~ gene) +
        labs(x = "Presence of mutation",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("0", "1")))

    p2 = ggplot(both %>% filter(mut == 1), aes(x=aneup_class, y=1/tot, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        guides(color="none") +
        facet_wrap(~ gene) +
        scale_color_manual(values=setNames(c("#cc9933", "#5500aa"), c("(0,0.1]", "(0.1,Inf]")), name="Aneuploidy") +
        labs(x="Aneuploidy", y="As fraction of total mutational load") +
        scale_y_log10() +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("(0,0.1]", "(0.1,Inf]")))

    p1 | p2
}

sys$run({
    aneup = get_aneup_purity()
    muts = get_muts(aneup)
    cnas = get_cnas(aneup)

    (plot_muts(muts) / plot_cnas(cnas)) + plot_layout(guides="collect")
})
