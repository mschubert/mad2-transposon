library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')

sc_wgs = function() {
    smps = c("401t", "419t", "413s")
    scs = file.path("../data/wgs", paste0(smps, ".rds")) %>%
        lapply(readRDS) %>%
        setNames(smps)
    plt$genome$heatmap_aneuHMM(scs) + guides(fill = guide_legend(title="Copy number"))
}

sys$run({
    args = sys$cmd$parse(
#        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
#        opt('p', 'poisson', 'rds', '../cis_analysis/poisson.rds'),
#        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
#        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
#        opt('r', 'rna_ins', 'txt', '../data/rnaseq_imfusion/insertions.txt'),
        opt('p', 'plotfile', 'pdf', 'FigS2-CIS.pdf')
    )

#    aneup = readRDS(args$aset)$meta %>%
#        select(genotype, sample, type, aneuploidy) %>%
#        mutate(type = factor(type))
#    levels(aneup$type)[levels(aneup$type) == "Other"] = "B-like"
#    ext = readRDS(args$ext)

    # insertion processing
#    rna_ins = readr::read_tsv(args$rna_ins, header=TRUE)
#    cis = readRDS(args$poisson)

    # create plot objects
    sc_wgs = sc_wgs()

    asm = sc_wgs + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 10, 6)
    print(asm)
    dev.off()
})
