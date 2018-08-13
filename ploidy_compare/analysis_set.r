library(dplyr)
library(cowplot)
library(patchwork)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
aneuf = import('tools/aneufinder')

#' source sample IDs are both 123S and T567
fix_ids = function(x) {
    letter = tolower(sub("[0-9]+", "", x))
    nums = sub("[ST]", "", x)
    paste0(nums, letter)
}

plot_comparison = function(aneups) {
    dens = ggplot(aneups, aes(x=aneuploidy)) +
        geom_density(aes(fill=type), alpha=0.5) +
        theme(axis.title.x=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank()) +
        guides(fill=FALSE)

    samples = ggplot(aneups, aes(x=aneuploidy, y=sample)) +
        geom_segment(aes(xend=aneuploidy, yend=sample), x=0,
                     color="lightgrey") +
        geom_point(aes(color=type, shape=tissue)) +
        theme(axis.text.y = element_text(size=8))

    dens + samples + plot_layout(ncol=1, heights=c(1,12))
}

plot_paircor = function(aneups) {
    paircor = aneups %>%
        group_by(sample, type, tissue) %>%
        summarize(aneuploidy = mean(aneuploidy)) %>%
        tidyr::spread("type", "aneuploidy") %>%
        GGally::ggpairs(columns=3:ncol(.), aes(shape=tissue))
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'metadata table', '../data/meta/meta.tsv'),
        opt('d', 'dna_seq', 'wgs ploidy', '../data/wgs/30cellseq.RData'),
        opt('r', 'rna_seq', 'eT ploidy', '../ploidy_from_rnaseq/eT_ploidy.RData'),
        opt('g', 'merge', 'fractions of high/low', 'analysis_set_merge.tsv'),
        opt('o', 'outfile', '.RData results', 'analysis_set.RData'),
        opt('p', 'plotfile', 'pdf', 'analysis_set.pdf'))

    meta = io$read_table(args$meta, header=TRUE)
    dna = io$load(args$dna_seq) %>%
        seq$aneuploidy(sample="Sample", assembly="GRCm38") %>%
        dplyr::rename(sample = Sample)
    rna = io$load(args$rna_seq)$segments %>%
        seq$aneuploidy(sample="sample", ploidy="expr", assembly="GRCm38") %>%
        mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))

    dna_merge = readr::read_tsv(args$merge) %>%
        left_join(dna %>% select(subset=sample, aneuploidy)) %>%
        group_by(sample) %>%
        summarize(aneuploidy = weighted.mean(aneuploidy, weight))
    dna_label = dna %>%
        mutate(sample = sub("(-high)|(-low)", "", sample),
               sample = paste0(sub("[^ST]+", "", sample), sub("[^0-9]+", "", sample)))

    sc = c("T401", "S413", "T419")
    sc_wgs = file.path("../data/wgs", paste0(sc, ".RData")) %>%
        io$load() %>%
        lapply(aneuf$consensus_ploidy) %>%
        dplyr::bind_rows(.id="sample") %>%
        seq$aneuploidy(sample="sample", width="length", assembly="GRCm38")

    tissues = setNames(c("spleen", "thymus"), c("S","T"))

    aneups = list(`WGS (merged)` = dna_merge,
                  `WGS (30-cell)` = dna_label,
                  `WGS (single-cell)` = sc_wgs,
                  `RNA-seq (eT)` = rna) %>%
        dplyr::bind_rows(.id="type") %>%
        filter(!is.na(aneuploidy)) %>%
        mutate(sample = forcats::fct_reorder(sample, aneuploidy),
               tissue = tissues[sub("Healthy|[0-9]+", "", sample)])

    merged = aneups %>%
        select(-coverage, -tissue) %>%
        filter(!duplicated(data.frame(type, sample)),
               !sample %in% c("S", "T")) %>%
        mutate(sample = fix_ids(sample)) %>%
        tidyr::spread("type", "aneuploidy") %>%
        group_by(sample) %>%
        summarize(aneup = `WGS (merged)` %or% `WGS (30-cell)` %or% `RNA-seq (eT)`) %>%
        arrange(-aneup) %>%
        inner_join(meta)

    pdf(9, 16, file=args$plotfile)
    print(plot_comparison(aneups))
    print(plot_paircor(aneups))
    dev.off()

    save(merged, file=args$outfile)
}
