library(dplyr)
library(ggplot2)
library(patchwork)
tcga = import('data/tcga')
sys = import('sys')

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

sys$run({
})
