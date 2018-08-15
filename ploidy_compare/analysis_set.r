library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')
seq = import('seq')
sys = import('sys')

plot_comparison = function(aneups, meta) {
    dset = aneups %>%
        inner_join(meta %>% select(sample, tissue, chosen=aneup_src)) %>%
        mutate(sample = forcats::fct_reorder(sample, aneuploidy),
               chosen = chosen == aneup_src)

    dens = ggplot(dset, aes(x=aneuploidy)) +
        geom_density(aes(fill=aneup_src), alpha=0.5) +
        theme(axis.title.x=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank()) +
        guides(fill=FALSE)

    samples = ggplot(dset, aes(x=aneuploidy, y=sample)) +
        geom_segment(aes(xend=aneuploidy, yend=sample), x=0, color="lightgrey") +
        geom_point(aes(color=aneup_src, shape=tissue, alpha=chosen), size=4) +
        scale_alpha_manual(values=c(0.4, 1)) +
        theme(axis.text.y = element_text(size=8))

    dens + samples + plot_layout(ncol=1, heights=c(1,12))
}

plot_paircor = function(aneups) {
    paircor = aneups %>%
        mutate(tissue = sub("[0-9]+", "", sample)) %>%
        group_by(sample, aneup_src, tissue) %>%
        summarize(aneuploidy = mean(aneuploidy)) %>%
        tidyr::spread("aneup_src", "aneuploidy") %>%
        GGally::ggpairs(columns=3:ncol(.), aes(shape=tissue))
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'metadata table', '../data/meta/meta.tsv'),
        opt('d', 'dna_seq', 'wgs ploidy', '../data/wgs/30cellseq.RData'),
        opt('r', 'rna_seq', 'eT ploidy', '../ploidy_from_rnaseq/eT_ploidy.RData'),
        opt('g', 'merge', 'fractions of high/low', 'analysis_set_merge.tsv'),
        opt('s', 'sc_seq', 'merged RData', '../data/wgs/sc_merge.RData'),
        opt('o', 'outfile', '.RData results', 'analysis_set.RData'),
        opt('p', 'plotfile', 'pdf', 'analysis_set.pdf'))

    meta_old = io$read_table(args$meta, header=TRUE)
    dna = io$load(args$dna_seq)$segments %>%
        seq$aneuploidy(sample="sample", assembly="GRCm38")
    rna = io$load(args$rna_seq)$segments %>%
        seq$aneuploidy(sample="sample", ploidy="ploidy", assembly="GRCm38")
    sc_wgs = io$load(args$sc_seq) %>%
        seq$aneuploidy(sample="sample", width="length", assembly="GRCm38")
    dna_merge = readr::read_tsv(args$merge) %>%
        left_join(dna %>% select(subset=sample, aneuploidy)) %>%
        group_by(sample) %>%
        summarize(aneuploidy = weighted.mean(aneuploidy, weight))

    aneups = list(`WGS (merged)` = dna_merge,
                  `WGS (30-cell)` = mutate(dna, sample=sub("-.*", "", sample)),
                  `WGS (single-cell)` = sc_wgs,
                  `RNA-seq (eT)` = rna) %>%
        dplyr::bind_rows(.id="aneup_src") %>%
        filter(!is.na(aneuploidy))

    # add HealthyS/T to meta?
    priority = c("WGS (merged)", "WGS (30-cell)", "RNA-seq (eT)") # no sc_wgs
    meta = aneups %>%
        mutate(ord = factor(aneup_src, levels=rev(priority), ordered=TRUE)) %>%
        group_by(sample) %>% top_n(1, ord) %>% ungroup() %>%
        select(-ord, -coverage) %>%
        inner_join(meta_old)

    pdf(8, 10, file=args$plotfile)
    print(plot_comparison(aneups, meta))
    print(plot_paircor(aneups))
    dev.off()

    save(meta, file=args$outfile)
})
