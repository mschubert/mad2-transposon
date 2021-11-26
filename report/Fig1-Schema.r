library(survival)
library(survminer)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
seq = import('seq')

cohort = function() {
    schema = grid::rasterGrob(magick::image_read("external/MouseCohort.svg"))
    ggplot() + annotation_custom(schema) + theme(panel.background=element_blank())
}

surv1 = function(meta_all) {
    gt_col = c("Mad2"="#e2e2e2", "Mad2 PB UT"="#d6e9ee", "Mad2 PB"="#afdde9", "PB"="#ffccaa")
    gt_col[] = shades::lightness(gt_col, shades::delta(-10))
    gt_col[] = shades::saturation(gt_col, shades::delta(0.2))
    gt_col[1] = "#b0b0b0"
    meta2 = meta_all %>%
        mutate(gtpc = factor(paste(genotype, pipc)),
               status = 1)
    levels(meta2$gtpc) = names(gt_col)
    fit = survfit(Surv(months_death, status) ~ gtpc, data=meta2)
    p1 = ggsurvplot(fit, data=meta2, palette=unname(gt_col))$plot +
        xlab("Months") +
        ylab("OS") +
        theme(legend.position="none")

    tdf = meta2 %>%
        group_by(gtpc) %>%
        summarize(n = n(), frac = n() / nrow(.)) %>%
        mutate(csum = rev(cumsum(rev(frac))),
               pos = frac/2 + lead(csum, 1),
               pos = if_else(is.na(pos), frac/2, pos))
    p2 = ggplot(tdf, aes(x=1, y=frac, fill=gtpc)) +
        geom_bar(stat="identity", width=1, color="white") +
        scale_y_continuous(breaks=tdf$pos, labels=tdf$n) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=gt_col, name="Genotype") +
        theme_void() +
        theme(axis.text.x = element_text()) +
        plot_layout(tag_level="new")

    (p1 | p2) + plot_layout(widths=c(2,1))
}

surv2 = function(meta) {
    meta = meta %>% mutate(status = 1)
    fit = survfit(Surv(months_death, status) ~ type, data=meta)
    p1 = ggsurvplot(fit, data=meta)$plot +
        xlab("Months") +
        ylab("LFS") +
        theme(legend.position="none")

    tdf = meta %>%
        group_by(type) %>%
        summarize(n = n(), frac = n() / nrow(.)) %>%
        mutate(csum = rev(cumsum(rev(frac))),
               pos = frac/2 + lead(csum, 1),
               pos = if_else(is.na(pos), frac/2, pos))
    p2 = ggplot(tdf, aes(x=1, y=frac, fill=type)) +
        geom_bar(stat="identity", width=1, color="white") +
        scale_y_continuous(breaks=tdf$pos, labels=tdf$n) +
        coord_polar("y", start=0) +
        guides(fill = guide_legend(title="Cancer type")) +
        theme_void() +
        theme(axis.text.x = element_text()) +
        plot_layout(tag_level="new")

    (p1 | p2) + plot_layout(widths=c(2,1))
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
        geom_segment(data=gen, aes(x=1, xend=len, y=0, yend=0), size=2.2, alpha=0.4, lineend="round") +
        geom_point(aes(x=start, y=0)) +
        geom_text(aes(x=start, label=external_gene_name), y=0, vjust=2, size=4) +
        facet_grid(. ~ chromosome_name, scales="free", space="free") +
        coord_cartesian(clip="off") +#, expand=FALSE) +
        theme_void() +
        theme(strip.background = element_blank(),
              strip.text.x = element_text(size=10),
              panel.spacing.x = unit(1, "mm"))
}

chroms = function(segs, meta) {
    wgs30 = segs %>%
        mutate(seqnames = factor(seqnames, levels=c(1:19,'X')),
               sample = factor(sample, levels=rev(levels(meta$sample)))) %>%
        inner_join(meta %>% select(sample, aneuploidy)) %>%
        mutate(copy.number = factor(round(ploidy)))

    plt$genome$heatmap(wgs30, cell="sample", y_dodge=2) +
        guides(fill=guide_legend(title="Copy number")) +
        theme(plot.background = element_rect(fill="transparent", color=NA),
              panel.background = element_rect(fill="transparent", color=NA),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              panel.spacing.x = unit(1, "mm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1),
              axis.text.y = element_text(size=8)) +
        coord_cartesian(clip="off") +
        plot_layout(tag_level="new")
}

genotype_weights = function(meta) {
    meta = meta %>%
        mutate(sw = factor(switching > 70))
    tw = meta %>%
        select(sample, spleen_g, thymus_g) %>%
        tidyr::gather("tissue", "weight", -sample) %>%
        mutate(tissue = sub("_g$", "", tissue))
    ct = ggplot(meta, aes(y=sample)) +
        geom_tile(aes(x="type", fill=type), color="black", size=0.2) +
#        scale_fill_manual(values=c("f"="#7570b3", "m"="#e6ab02")) +
        coord_fixed(clip="off", expand=FALSE) +
        guides(fill=FALSE)
    sex = ggplot(meta %>% mutate(Sex=sex, sex=factor("sex")), aes(y=sample)) +
        geom_tile(aes(x=sex, fill=Sex), color="black", size=0.2) +
        scale_fill_manual(values=c("f"="#7570b3", "m"="#e6ab02")) +
        coord_fixed(clip="off", expand=FALSE)
    gt = ggplot(meta %>% mutate(gt=genotype, genotype=factor("genotype")), aes(y=sample)) +
        geom_tile(aes(x=genotype, fill=gt), color="black", size=0.2) +
        guides(fill=guide_legend(title="Genotype")) +
        scale_fill_manual(values=c("Mad2 PB Mx1-Cre"="darkorchid", "PB Mx1-Cre"="white")) +
        coord_fixed(clip="off", expand=FALSE)
    sw = ggplot(meta, aes(x=sample, y=pmax(switching, 10))) +
        geom_col(aes(fill=sw), color="white") +
        scale_fill_manual(values=c("TRUE"="purple", "FALSE"="grey"), guide=FALSE) +
        coord_flip(expand=FALSE, clip="off") +
        scale_y_continuous(breaks=50, labels="50% switch") +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank())
    gtsw = (gt | sw) + plot_layout(widths=c(1,2))
    tumors = ggplot(tw, aes(x=tissue, y=sample)) +
        geom_point(aes(size=weight), alpha=0.2) +
        coord_cartesian(clip="off") +
        guides(size=guide_legend(title="Weight (grams)")) +
        scale_size_area(breaks=c(0.5,1)) &
        theme_minimal()
    ct + sex + gtsw + tumors + plot_layout(widths=c(1,1,3,3), tag_level="new") &
        theme(plot.margin = margin(0,0,0,0,"mm"),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(angle=60, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'rds', '../data/meta/meta.rds'),
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig1-Schema.pdf')
    )

    meta_all = readRDS(args$meta)

    aset = readRDS(args$aset)
    meta = aset$meta %>%
        filter(!is.na(aneuploidy)) %>%
        mutate(sample = forcats::fct_reorder(sample, aneuploidy))
    segs = aset$segs %>%
        filter(is_pref_src)

    cohort_ov = cohort() & theme(plot.margin = margin(0,-5,0,0,"cm"))
    topright = wrap_plots(wrap_elements((surv1(meta_all) / surv2(meta)) & theme(plot.margin = margin(0,0,-2.5,-4,"mm"))))
    top = (cohort_ov | topright) + plot_layout(widths=c(1.6,1))
    btm = chrom_genes() + plot_spacer() + chroms(segs, meta) + genotype_weights(meta) +
        plot_layout(widths=unit(c(1,3), c("null","cm")), heights=c(1,18), guides="collect")

    asm = wrap_plots(top) + wrap_plots(wrap_elements(btm)) +
        plot_annotation(tag_levels='a') + plot_layout(ncol=1, heights=c(1.2,2)) &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 15, 12)
    print(asm)
    dev.off()
})
