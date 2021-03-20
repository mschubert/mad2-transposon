library(dplyr)
library(ggplot2)
library(patchwork)
library(plyranges)
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
    merge_interval = function(smp) {
        if (length(unique(smp$subset)) == 1)
            return(smp %>% select(seqnames, start, end, copy.number))
        smp = GenomicRanges::makeGRangesListFromDataFrame(smp, "subset", keep.extra.columns=TRUE)
        plyranges::join_overlap_intersect(smp[[1]], smp[[2]]) %>%
            as.data.frame() %>% as_tibble() %>%
            mutate(copy.number = (copy.number.x * weight.x + copy.number.y * weight.y) / (weight.x + weight.y)) %>%
            select(seqnames, start, end, copy.number)
    }
    wgs30 = wgs$segments %>%
        as_tibble() %>%
        filter(seqnames %in% c(1:19,'X')) %>%
        dplyr::rename(subset = sample) %>%
        left_join(wgs_merge %>% select(-comment), by="subset") %>%
        mutate(sample = ifelse(is.na(sample), subset, sample)) %>%
        group_by(sample) %>%
            tidyr::nest() %>%
            mutate(joined = purrr::map(data, merge_interval)) %>%
        ungroup() %>%
        select(-data) %>%
        tidyr::unnest("joined") %>%
        inner_join(aset %>% select(sample, aneuploidy)) %>%
        mutate(copy.number = factor(round(copy.number)),
               sample = forcats::fct_reorder(sample, aneuploidy))

    ggplot(wgs30, aes(y=sample, yend=sample)) +
        geom_segment(aes(x=start, xend=end, color=copy.number), size=1.5) +
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
    gt + tumors + plot_layout(widths=c(1,2), tag_level="new") &
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

    asm = ((cohort() | surv(meta) | types(meta)) + plot_layout(widths=c(1.8,1,1))) /
#        plt$text("pathology imgs go here") +
        (chrom_genes() + plot_layout(widths=c(5,1)) + plot_spacer() +
         chroms(wgs, aset, wgs_merge) + genotype_weights(meta) +
            plot_layout(widths=c(10,1), heights=c(1,50), guides="collect")) +
        plot_annotation(tag_levels='a') + plot_layout(heights=c(1,2)) &
        theme(plot.margin=margin(0.25, 0.25, 0.25, 0.25, "mm"),
              plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 15, 10)
    print(asm)
    dev.off()
})
