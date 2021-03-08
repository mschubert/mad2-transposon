library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')

cohort = function() {
    schema = grid::rasterGrob(magick::image_read("external/MouseCohort.svg"))
    splot1 = ggplot() + annotation_custom(schema)
    schema2 = grid::rasterGrob(magick::image_read("external/MouseCohort_process.svg"))
    splot2 = ggplot() + annotation_custom(schema2)
    schema3 = grid::rasterGrob(magick::image_read("external/MouseCohort_seq.svg"))
    splot3 = ggplot() + annotation_custom(schema3)

    splot1 + splot2 + splot3 + plot_layout(ncol=1, heights=c(3,0.8,1.1)) +
        plot_annotation(tag_levels='a')
}

surv = function() {
}

chroms = function() {
    mrg = readr::read_tsv("../ploidy_compare/analysis_set_merge.tsv")
    wgs30 = readRDS(args$wgs)$segments
    aneup = seq$aneuploidy(wgs30, assembly="GRCm38", sample="sample") %>%
        arrange(aneuploidy) %>%
        filter(sample %in% c(aset$sample, mrg$subset)) %>%
        mutate(sample = factor(sample, levels=sample))
}

genotype_weights = function(meta) {
    tw = meta %>%
        select(sample, spleen_g, thymus_g) %>%
        tidyr::gather("tissue", "weight", -sample) %>%
        mutate(tissue = sub("_g$", "", tissue))
    gt = ggplot(meta %>% mutate(gt=genotype, genotype=factor("genotype")), aes(y=sample)) +
        geom_tile(aes(x=genotype, fill=gt))
    tumors = ggplot(tw, aes(x=tissue, y=sample)) +
        geom_point(aes(size=weight), alpha=0.7) +
        scale_size_area()
    gt + tumors + guide_area() + plot_layout(widths=c(1,2,8), guides="collect") &
        theme_minimal() &
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'tsv', '../data/meta/meta.tsv'),
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('w', 'wgs', 'rds', '../data/wgs/30cellseq.rds'),
#        opt('', '', '', ''),
        opt('p', 'plotfile', 'pdf', 'Fig1-Schema.pdf')
    )


    meta = readr::read_tsv(args$meta)

})


models = models[levels(aneups$sample)]
for (i in seq_along(models)) {
    int = as.integer(round(gr[[models[[i]]$ID]]$ploidy))
    cn = paste0(int, "-somy")
    cn = factor(cn, levels=levels(models[[i]]$segments$state))
    models[[i]]$segments$state = cn
    models[[i]]$segments$copy.number = int
}
