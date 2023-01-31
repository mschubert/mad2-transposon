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

get_muts = function(aneup) {
    nosyns = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
        "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site") # "Splice_Region",

    load_fl = function(coh) tcga$mutations(coh) %>%
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol, vclass=Variant_Classification)
    m = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows() %>%
        group_by(cohort, Sample) %>%
            mutate(n = n_distinct(gene)) %>%
            summarize(tot = unique(n),
                      STAT1 = as.integer("STAT1" %in% gene[vclass %in% nosyns]),
                      TP53 = as.integer("TP53" %in% gene[vclass %in% nosyns]),
                      JAK1 = as.integer("JAK1" %in% gene[vclass %in% nosyns]),
#                      B2M = as.integer("B2M" %in% gene[vclass %in% nosyns]),
                      TTN = as.integer("TTN" %in% gene[vclass %in% nosyns])) %>%
        ungroup()
    aneup %>%
        inner_join(m) %>%
        tidyr::gather("gene", "mut", -(Sample:tot)) %>%
        filter(tot <= 5000)
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
                      `IFNA1 loss` = as.integer(gistic[gene == "IFNA1"] == -2),
                      `IFNB1 loss` = as.integer(gistic[gene == "IFNB1"] == -2),
                      `CDKN2A loss` = as.integer(gistic[gene == "CDKN2A"] == -2)) %>%
        ungroup()
    inner_join(aneup, cfrac) %>%
        tidyr::gather("gene", "cna", -(Sample:n_del)) #%>%
#        mutate(gene = factor(gene, levels=c("MYC amp", "IFNA1 loss", "TP53 loss")))
}

plot_muts = function(muts) {
    ev_fill = c(Absent="#b3cde3", Present="#fbb4ae")
    cna_cols = c(`(0,0.1]`="#a65628", `(0.1,Inf]`="#5500aa")
    muts$event = ifelse(muts$mut, "Present", "Absent") %>% factor(levels=names(ev_fill))

    p1 = ggplot(muts, aes(x=event, y=aneup, fill=event)) +
        geom_boxplot(outlier.shape=NA) +
        scale_fill_manual(values=ev_fill, name="Mutation") +
        facet_wrap(~ gene) +
        labs(x = "Presence of mutation",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(names(ev_fill)))

    #todo: fix y_position, tip_length (ggsignif/issues/48)
    p2 = ggplot(muts, aes(x=aneup_class, y=tot, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA, fill="#dedede") +
        coord_cartesian(ylim=c(NA, 450)) +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        labs(x = "Aneuploidy",
             y = "Total number of mutated genes") +
        ggsignif::geom_signif(color="black", test=wilcox.test, tip_length=1e-3,
            comparisons=list(names(cna_cols)), y_position=200)

    p3 = ggplot(muts %>% filter(mut == "1"), aes(x=aneup_class, y=1/tot, color=aneup_class)) +
        geom_boxplot(aes(fill=event), outlier.shape=NA) +
        scale_color_manual(values=cna_cols, guide="none") +
        scale_fill_manual(values=ev_fill, guide="none") +
        facet_wrap(~ gene) +
        labs(x="Aneuploidy", y="As fraction of mutated genes") +
        scale_y_log10() +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(names(cna_cols)))

    (p1 | p2 | p3) + plot_layout(widths=c(2,1,2))
}

plot_cnas = function(cnas) {
    ev_fill = c(Absent="#ccebc5", Present="#fed9a6")
    cna_cols = c(`(0,0.1]`="#a65628", `(0.1,Inf]`="#5500aa")
    cnas$event = ifelse(cnas$cna, "Present", "Absent") %>% factor(levels=names(ev_fill))

    p1 = ggplot(cnas, aes(x=event, y=aneup, fill=event)) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_fill_manual(values=ev_fill, name="Deletion") +
        labs(x = "Presence of event",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(names(ev_fill)))

    #todo: fix y_position, tip_length (ggsignif/issues/48)
    p2 = ggplot(cnas, aes(x=aneup_class, y=n_del, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA, fill="#dedede") +
        coord_cartesian(ylim=c(NA, 260)) +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        labs(x = "Aneuploidy",
             y = "Total number of lost genes") +
        ggsignif::geom_signif(color="black", test=wilcox.test, tip_length=1e-3,
            comparisons=list(names(cna_cols)), y_position=10)

    p3 = ggplot(cnas, aes(x=aneup_class, y=cna/(n_del), color=aneup_class)) +
        geom_boxplot(aes(fill=event), outlier.shape=NA) +
        scale_color_manual(values=cna_cols, guide="none") +
        scale_fill_manual(values=ev_fill, guide="none") +
        facet_wrap(~ gene) +
        scale_y_log10() +
        labs(x = "Aneuploidy",
             y = "As fraction of genes lost") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(names(cna_cols)))

    (p1 | p2 | p3) + plot_layout(widths=c(2,1,2))
}

sys$run({
    aneup = get_aneup_purity()
    muts = get_muts(aneup)
    cnas = get_cnas(aneup)

    asm = (plot_muts(muts) / plot_cnas(cnas)) + plot_layout(guides="collect") &
        theme_classic() & theme(
            strip.background = element_blank(),
            strip.text.x = element_text(size=12, face="bold")
        )

    pdf("mut_frac.pdf", 12, 8)
    print(asm)
    dev.off()

    saveRDS(list(aneup=aneup, muts=muts, cnas=cnas), file="mut_frac.rds")
})
