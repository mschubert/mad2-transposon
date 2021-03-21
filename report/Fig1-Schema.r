library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
seq = import('seq')

cohort = function() {
    schema = grid::rasterGrob(magick::image_read("external/MouseCohort.svg"))
    ggplot() + annotation_custom(schema)
}

surv = function(meta) {
    surv = grid::rasterGrob(magick::image_read("external/Overall survival in months - age.pdf"))
    ggplot() + annotation_custom(surv) + theme_void()
}

types = function(meta) {
    tdf = meta %>%
        group_by(type) %>%
        summarize(frac = n() / nrow(.))
    types = ggplot(tdf, aes(x="", y=frac, fill=type)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
#        geom_text(aes(y=ypos, label=frac), color = "white", size=6) +
        scale_fill_brewer(palette="Set1") +
        theme_void()
}

chrom_genes = function() {
    gen = seq$genome(assembly="GRCm38") %>%
        GenomeInfoDb::seqlengths()
    gen = data.frame(chromosome_name=c(1:19,'X'), len=gen[c(1:19,'X')]/1e6) %>%
        mutate(chromosome_name = factor(chromosome_name, levels=c(1:19,'X')))
    hlg = seq$coords$gene(dset="mmusculus_gene_ensembl", assembly="GRCm38") %>%
        filter(external_gene_name %in% c("Trp53", "Ets1", "Erg", "Pten", "Notch1", "Myc", "Kras")) %>%
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

chroms = function(segs, meta) {
    wgs30 = segs %>%
        mutate(seqnames = factor(seqnames, levels=c(1:19,'X')),
               sample = factor(sample, levels=levels(meta$sample))) %>%
        inner_join(meta %>% select(sample, aneuploidy)) %>%
        mutate(copy.number = factor(round(ploidy)),
               cell = sample)

    plt$genome$heatmap(wgs30) +
        guides(fill=guide_legend(title="Copy number")) +
        theme(plot.background = element_rect(fill="transparent", color=NA),
              panel.background = element_rect(fill="transparent", color=NA),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              panel.spacing.x = unit(0.5, "mm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1),
              axis.text.y = element_text(size=5.5)) +
        coord_cartesian(clip="off", expand=FALSE)
}

genotype_weights = function(meta) {
    tw = meta %>%
        select(sample, spleen_g, thymus_g) %>%
        tidyr::gather("tissue", "weight", -sample) %>%
        mutate(tissue = sub("_g$", "", tissue))
    sex = ggplot(meta %>% mutate(Sex=sex, sex=factor("sex")), aes(y=sample)) +
        geom_tile(aes(x=sex, fill=Sex), color="black", size=0.2) +
        scale_fill_manual(values=c("f"="#7570b3", "m"="#e6ab02")) +
        coord_fixed(clip="off")
    gt = ggplot(meta %>% mutate(gt=genotype, genotype=factor("genotype")), aes(y=sample)) +
        geom_tile(aes(x=genotype, fill=gt), color="black", size=0.2) +
        guides(fill=guide_legend(title="Genotype")) +
        scale_fill_manual(values=c("Mad2 PB Mx1-Cre"="darkorchid", "PB Mx1-Cre"="white")) +
        coord_fixed(clip="off")
    tumors = ggplot(tw, aes(x=tissue, y=sample)) +
        geom_point(aes(size=weight), alpha=0.2) +
        coord_cartesian(clip="off") +
        guides(size=guide_legend(title="Weight (grams)")) +
        scale_size_area()
    #todo: mad2 switching in % as bar?
    sex + gt + tumors + plot_layout(widths=c(1,1,3), tag_level="new") &
        theme_minimal() &
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
#        opt('', '', '', ''),
        opt('p', 'plotfile', 'pdf', 'Fig1-Schema.pdf')
    )

    aset = readRDS(args$aset)
    meta = aset$meta %>%
        filter(!is.na(aneuploidy)) %>%
        mutate(sample = forcats::fct_reorder(sample, aneuploidy))
    segs = aset$segs %>%
        filter(is_pref_src)

    asm = ((cohort() | surv(meta) | types(meta)) + plot_layout(widths=c(1.8,1,1))) /
#        plt$text("pathology imgs go here") +
        (chrom_genes() + plot_layout(widths=c(5,1)) + plot_spacer() +
         chroms(segs, meta) + genotype_weights(meta) +
            plot_layout(widths=c(10,1), heights=c(1,50), guides="collect")) +
        plot_annotation(tag_levels='a') + plot_layout(heights=c(1,2)) &
        theme(plot.margin=margin(0.25, 0.25, 0.25, 0.25, "mm"),
              plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 15, 10)
    print(asm)
    dev.off()
})
