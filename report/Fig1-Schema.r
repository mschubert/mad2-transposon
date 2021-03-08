library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
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

surv = function(meta) {
    surv = grid::rasterGrob(magick::image_read("external/Overall survival in months - age.pdf"))
    splot = ggplot() + annotation_custom(surv)

    tdf = meta %>%
        group_by(type) %>%
        summarize(frac = n() / nrow(.))
    types = ggplot(tdf, aes(x="", y=frac, fill=type)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
#        geom_text(aes(y=ypos, label=frac), color = "white", size=6) +
        scale_fill_brewer(palette="Set1")

    splot / types & theme_void()
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
        geom_point(aes(x=genotype, fill=gt), size=2, shape=22)
    tumors = ggplot(tw, aes(x=tissue, y=sample)) +
        geom_point(aes(size=weight), alpha=0.7) +
        coord_cartesian(clip="off") +
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

    pdf(args$plotfile, 16, 19)
    ((cohort() | surv(meta)) + plot_layout(widths=c(3,2))) /
#        plt$text("x goes here") +
        genotype_weights(meta) +
        plot_annotation(tag_levels='a') +
        plot_layout(heights=c(1,2))
    dev.off()
})
