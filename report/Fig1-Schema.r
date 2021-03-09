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

chrom_genes = function() {
    gen = seq$genome(assembly="GRCm38") %>%
        GenomeInfoDb::seqlengths()
    gen = data.frame(chromosome_name=c(1:19,'X'), len=gen[c(1:19,'X')]/1e6) %>%
        mutate(chromosome_name = factor(chromosome_name, levels=c(1:19,'X')))
    hlg = seq$coords$gene(dset="mmusculus_gene_ensembl", assembly="GRCm38") %>%
        filter(external_gene_name %in% c("Trp53", "Ets1", "Erg", "Pten", "Notch1", "Myc")) %>%
        mutate(start = start_position / 1e6,
               end = end_position / 1e6,
               chromosome_name = factor(chromosome_name, levels=c(1:19,'X')))
    ggplot(hlg) +
        geom_segment(data=gen, aes(x=1, xend=len, y=0, yend=0), size=2, alpha=0.4, lineend="round") +
        geom_point(aes(x=start, y=0)) +
        geom_text(aes(x=start, label=external_gene_name), y=0, vjust=2, size=3) +
        facet_grid(. ~ chromosome_name, scales="free", space="free") +
        coord_cartesian(clip="off") +#, expand=FALSE) +
        theme_void() +
        theme(strip.background = element_blank(),
              panel.spacing.x = unit(1, "mm"))
}

chroms = function(wgs, aset, wgs_merge) {
    wgs30 = wgs$segments %>%
        filter(seqnames %in% c(1:19,'X'),
               state %in% paste0(c(1:5), "-somy")) %>%
        mutate(start = start / 1e6,
               end = end / 1e6)

    aneup = seq$aneuploidy(wgs30, assembly="GRCm38", sample="sample") %>%
        arrange(aneuploidy) %>%
        filter(sample %in% c(aset$sample, wgs_merge$subset)) %>%
        mutate(sample = factor(sample, levels=sample))

    ggplot(wgs30, aes(y=sample, yend=sample)) +
        geom_segment(aes(x=start, xend=end, color=state), size=1.5) +
        facet_grid(. ~ seqnames, scales="free", space="free") +
        scale_x_continuous(breaks=c(50, 100, 150)) +
        theme(plot.background = element_rect(fill="transparent", color=NA),
              panel.background = element_rect(fill="transparent", color=NA),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              panel.spacing.x = unit(0.5, "mm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1)) +
        coord_cartesian(clip="off", expand=FALSE) +
        xlab("Position (Mb)")
}

genotype_weights = function(meta) {
    tw = meta %>%
        select(sample, spleen_g, thymus_g) %>%
        tidyr::gather("tissue", "weight", -sample) %>%
        mutate(tissue = sub("_g$", "", tissue))
    gt = ggplot(meta %>% mutate(gt=genotype, genotype=factor("genotype")), aes(y=sample)) +
        geom_point(aes(x=genotype, fill=gt), size=2.5, shape=22) +
        coord_cartesian(clip="off")
    tumors = ggplot(tw, aes(x=tissue, y=sample)) +
        geom_point(aes(size=weight), alpha=0.4) +
        coord_cartesian(clip="off") +
        scale_size_area()
    #todo: mad2 switching in % as bar?
    gt + tumors + plot_layout(widths=c(1,2)) &
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
        opt('m', 'wgs_merge', 'rds', '../ploidy_compare/analysis_set_merge.tsv'),
#        opt('', '', '', ''),
        opt('p', 'plotfile', 'pdf', 'Fig1-Schema.pdf')
    )

    meta = readr::read_tsv(args$meta)
    aset = readRDS(args$aset)
    wgs_merge = readr::read_tsv(args$wgs_merge)
    wgs = readRDS(args$wgs)

    pdf(args$plotfile, 15, 14)
    ((cohort() | surv(meta)) + plot_layout(widths=c(7,4))) /
#        plt$text("x goes here") +
        (chrom_genes() + plot_layout(widths=c(5,1)) + plot_spacer() +
         chroms(wgs, aset, wgs_merge) + genotype_weights(meta) +
            plot_layout(widths=c(5,1), heights=c(1,50), guides="collect")) +
        plot_annotation(tag_levels='a') + plot_layout(heights=c(2,1.5))
    dev.off()
})
