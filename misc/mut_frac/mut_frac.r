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
                      `IFNA1 loss` = as.integer(gistic[gene == "IFNA1"] == -2),
                      `IFNB1 loss` = as.integer(gistic[gene == "IFNB1"] == -2),
                      `CDKN2A loss` = as.integer(gistic[gene == "CDKN2A"] == -2)) %>%
        ungroup()
    inner_join(aneup, cfrac) %>%
        tidyr::gather("gene", "cna", -(Sample:n_del)) #%>%
#        mutate(gene = factor(gene, levels=c("MYC amp", "IFNA1 loss", "TP53 loss")))
}

plot_cnas = function(cnas) {
    cna_cols = setNames(c("#cc9933", "#5500aa"), c("(0,0.1]", "(0.1,Inf]"))

    p1 = ggplot(cnas, aes(x=factor(cna), y=aneup, color=factor(cna))) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_color_brewer(palette="Set1", name="Genomic\nevent") +
        labs(x = "Presence of event",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("0", "1")))

    #todo: fix y_position, tip_length (ggsignif/issues/48)
    p2 = ggplot(cnas, aes(x=aneup_class, y=n_del, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        coord_cartesian(ylim=c(NA, 260)) +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        labs(x = "Aneuploidy",
             y = "Total number of lost genes") +
        ggsignif::geom_signif(color="black", test=wilcox.test, tip_length=1e-3,
            comparisons=list(c("(0,0.1]", "(0.1,Inf]")), y_position=10)

    p3 = ggplot(cnas, aes(x=aneup_class, y=cna/(n_del), color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_y_log10() +
        labs(x = "Aneuploidy",
             y = "As fraction of genes lost") +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test,
            comparisons=list(c("(0,0.1]", "(0.1,Inf]")))

    (p1 | p2 | p3) + plot_layout(widths=c(2,1,2))
}

all_muts_volc = function(aneup) {
    plt = import('plot')
    gset = import('genesets')

    ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
    is_prot = GeneBiotypeFilter("protein_coding")
    is_chr = SeqNameFilter(c(1:22,'X','Y'))
    tx = transcripts(ens106, filter=c(is_prot, is_chr)) %>%
        plyranges::filter(tx_is_canonical == 1)
    ex = exonsBy(ens106)[names(tx)]
    names(ex) = sub("-[0-9]+$", "", tx$tx_external_name)
    ex = ex[!is.na(names(ex))] %>%
        lapply(. %>% summarize(glen = sum(end(.) - start(.))) %>% as.data.frame()) %>%
        bind_rows(.id="gene")
    filter = dplyr::filter

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
    f1 = gset$test_fet(res$gene[res$estimate>quantile(res$estimate,0.8)], sets[[1]], valid=res$gene)
    plt$volcano(f1 %>% mutate(estimate=log2(estimate)))
    s1 = gset$test_lm(res, stat="estimate", sets[[1]])
    plt$volcano(s1)
    f2 = gset$test_fet(res$gene[res$estimate>quantile(res$estimate,0.8)], sets[[2]], valid=res$gene, min=4)
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

dn_ds_try = function(aneup) {
    load_fl = function(coh) tcga$mutations(coh) %>%
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol, vclass=Variant_Classification)
    m2 = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows() %>%
        inner_join(aneup)

    nosyns = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
        "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site") # "Splice_Region",
    syns = c("Silent") #, "Intron") #"3'Flank", "3'UTR", "5'UTR", "5'Flank", "RNA", "IGR") #all: EP300 hit
    fracs = m2 %>%
        group_by(Sample) %>%
            filter(!is.na(aneup_class)) %>%
        group_by(aneup_class, gene) %>%
            summarize(dn = length(unique(Sample[vclass %in% nosyns])),
                      ds = length(unique(Sample[vclass %in% syns])))
    res2 = fracs %>%
        group_by(gene) %>%
            filter(n() == 2) %>%
            arrange(aneup_class) %>%
            summarize(res = list(broom::tidy(fisher.test(matrix(c(ds, dn), ncol=2))))) %>%
        tidyr::unnest(res) %>%
        dplyr::select(-method, -alternative) %>%
        arrange(p.value)
    # JAK1 evidence (p=0.02); STAT1, TP53 unchanged (p=1, 0.2)
    # with only Silent, no Intron: p=0.02; 1, 0.6
    # -> across IFNA all genes ~0.5 enrich (so opposite effect)
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

plot_muts = function(muts) {
    cna_cols = setNames(c("#cc9933", "#5500aa"), c("(0,0.1]", "(0.1,Inf]"))

    p1 = ggplot(muts, aes(x=factor(mut), y=aneup, color=factor(mut))) +
        geom_boxplot(outlier.shape=NA) +
        scale_color_brewer(palette="Set1", name="Genomic\nevent") +
        facet_wrap(~ gene) +
        labs(x = "Presence of mutation",
             y = "Aneuploidy") +
        ggsignif::geom_signif(color="black", test=wilcox.test, comparisons=list(c("0", "1")))

    #todo: fix y_position, tip_length (ggsignif/issues/48)
    p2 = ggplot(muts, aes(x=aneup_class, y=tot, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        coord_cartesian(ylim=c(NA, 450)) +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        labs(x = "Aneuploidy",
             y = "Total number of mutated genes") +
        ggsignif::geom_signif(color="black", test=wilcox.test, tip_length=1e-3,
            comparisons=list(c("(0,0.1]", "(0.1,Inf]")), y_position=200)

    p3 = ggplot(muts %>% filter(mut == 1), aes(x=aneup_class, y=1/tot, color=aneup_class)) +
        geom_boxplot(outlier.shape=NA) +
        facet_wrap(~ gene) +
        scale_color_manual(values=cna_cols, name="Aneuploidy") +
        labs(x="Aneuploidy", y="As fraction of mutated genes") +
        scale_y_log10() +
        ggsignif::geom_signif(color="black", test=wilcox.test,
            comparisons=list(c("(0,0.1]", "(0.1,Inf]")))

    (p1 | p2 | p3) + plot_layout(widths=c(2,1,2))
}

sys$run({
    aneup = get_aneup_purity()
    muts = get_muts(aneup)
    cnas = get_cnas(aneup)

    (plot_muts(muts) / plot_cnas(cnas)) + plot_layout(guides="collect") & theme_classic()
})
