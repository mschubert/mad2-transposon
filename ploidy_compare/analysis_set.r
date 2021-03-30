library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
seq = import('seq')
sys = import('sys')

#' Show ordered aneuploidies for all samples and sources
plot_comparison = function(aneups, meta) {
    smp_order = aneups %>% filter(is_pref_src) %>% arrange(aneuploidy) %>% pull(sample)
    dset = aneups %>%
        mutate(sample = factor(sample, levels=smp_order))

    dens = ggplot(dset, aes(x=aneuploidy)) +
        geom_density(aes(fill=ploidy_src), alpha=0.5) +
        theme(axis.title.x=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank()) +
        guides(fill=FALSE)

    samples = ggplot(dset, aes(x=aneuploidy, y=sample)) +
        geom_segment(aes(xend=aneuploidy, yend=sample), x=0, color="lightgrey") +
        geom_point(aes(color=ploidy_src, alpha=is_pref_src), size=4) + # shape=tissue
        scale_alpha_manual(values=c(0.3, 0.8)) +
        theme(axis.text.y = element_text(size=8))

    dens + samples + plot_layout(ncol=1, heights=c(1,12))
}

#' Pariwise correlation of aneuploidies from RNA-seq or 30-cell WGS
plot_paircor = function(aneups) {
    paircor = aneups %>%
        filter(ploidy_src %in% c("RNA-seq (eT)", "WGS (30-cell)")) %>%
        mutate(tissue = sub("[0-9]+", "", sample)) %>%
        group_by(sample, ploidy_src, tissue) %>%
        summarize(aneuploidy = mean(aneuploidy)) %>%
        tidyr::spread("ploidy_src", "aneuploidy") %>%
        GGally::ggpairs(columns=3:ncol(.), aes(shape=tissue))
}

#' Merge segments where we 30-cell DNA sequenced high and low populations
#'
#' @param segs       A data.frame of segments: sample, subset, seqnames, start, stop, ploidy
#' @param wgs_merge  Merge table from analysis_set_merge.tsv: sample, subset, weight
merge_pops = function(segs, wgs_merge) {
    merge_interval = function(smp) {
        if (length(unique(smp$subset)) == 1) # merge only 2 pops
            return(smp %>% select(seqnames, start, end, ploidy))
        smp = GenomicRanges::makeGRangesListFromDataFrame(smp, "subset", keep.extra.columns=TRUE)
        plyranges::join_overlap_intersect(smp[[1]], smp[[2]]) %>%
            as.data.frame() %>% as_tibble() %>%
            mutate(ploidy = (ploidy.x * weight.x + ploidy.y * weight.y) / (weight.x + weight.y)) %>%
            select(seqnames, start, end, ploidy)
    }

    as_tibble(segs) %>%
        left_join(wgs_merge, by=c("sample", "subset")) %>%
        group_by(sample) %>%
            tidyr::nest() %>%
            mutate(joined = purrr::map(data, merge_interval)) %>%
        ungroup() %>%
        select(-data) %>%
        tidyr::unnest("joined")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'metadata table', '../data/meta/meta.rds'),
        opt('d', 'dna_seq', 'wgs ploidy', '../data/wgs/30cellseq.rds'),
        opt('r', 'rna_seq', 'eT ploidy', '../ploidy_from_rnaseq/eT_ploidy.rds'),
        opt('g', 'merge', 'fractions of high/low', 'analysis_set_merge.tsv'),
        opt('s', 'sc_seq', 'merged rds', '../data/wgs/sc_merge.rds'),
        opt('o', 'outfile', 'rds results', 'analysis_set.rds'),
        opt('p', 'plotfile', 'pdf', 'analysis_set.pdf')
    )

    meta_old = readRDS(args$meta)
    rna = readRDS(args$rna_seq)$segments
    sc_wgs = readRDS(args$sc_seq) %>%
        mutate(start=1, end=length) %>% # todo: actually get segments + merge here
        filter(seqnames %in% c(1:19,'X')) # JH584304.1 has end=NA, issues with aneup calc

    dna = readRDS(args$dna_seq)$segments %>%
        as_tibble() %>%
        mutate(subset = ifelse(grepl("-(high|low)$", sample), sample, NA),
               sample = sub("-(high|low)$", "", sample))
    dna_merge = readr::read_tsv(args$merge)
    dna_smp = merge_pops(dna, dna_merge)

    priority = c("WGS (30-cell)", "WGS (single-cell)", "RNA-seq (eT)", "WGS (high+low)")
    segs = list(`WGS (30-cell)` = dna_smp,
                `WGS (high+low)` = dna %>% filter(!is.na(subset)),
                `WGS (single-cell)` = sc_wgs,
                `RNA-seq (eT)` = rna) %>%
        dplyr::bind_rows(.id="ploidy_src") %>%
        select(sample, ploidy_src, seqnames, start, end, ploidy) %>%
        mutate(ord = factor(ploidy_src, levels=priority, ordered=TRUE)) %>%
        group_by(sample) %>%
            mutate(is_pref_src = ord == min(ord)) %>%
        ungroup() %>%
        select(-ord)

    aneups = segs %>% mutate(smp_src = paste(sample, ploidy_src), width = end - start)
    nums = seq$aneuploidy(aneups, sample="smp_src", ploidy="ploidy", assembly="GRCm38")
    aneups = aneups %>%
        select(smp_src, sample, ploidy_src, is_pref_src) %>%
        distinct() %>%
        inner_join(nums, by="smp_src") %>%
        select(-smp_src)

    meta = meta_old %>%
        inner_join(aneups %>% filter(is_pref_src) %>% select(sample, aneuploidy, ploidy_src)) %>%
        mutate(`T-cell` = (type == "T-cell" & !is.na(type)) + 0,
               Myeloid = (type == "Myeloid" & !is.na(type)) + 0,
               `B-like` = (type == "B-like" & !is.na(type)) + 0) %>%
        select(sample, everything())

    pdf(8, 10, file=args$plotfile)
    print(plot_comparison(aneups))
    print(plot_paircor(aneups))
    dev.off()

    saveRDS(list(meta=meta, aneups=aneups, segs=segs), file=args$outfile)
})
