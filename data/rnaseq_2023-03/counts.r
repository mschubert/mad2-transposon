library(dplyr)
library(ggplot2)
sys = import('sys')

#' Read a character vector of STAR read count files
read_files = function(fnames) {
    samples = sub(".ReadsPerGene.out.tab", "", basename(fnames), fixed=TRUE)
    lapply(fnames, readr::read_tsv, col_names=FALSE) %>%
        setNames(samples) %>%
        bind_rows(.id="sample") %>%
        dplyr::rename(feature=X1, unstranded=X2, sense=X3, antisense=X4)
}

#' Summarize STAR read count statistics for different ways of counting
read_stats = function(reads) {
    reads %>%
        mutate(feature = ifelse(grepl("^N_", feature), feature, "N_mapped")) %>%
        group_by(sample, feature) %>%
            summarize_if(is.numeric, sum) %>%
        ungroup() %>%
        tidyr::pivot_longer(c(unstranded, sense, antisense), names_to="count_type", values_to="reads") %>%
        mutate(count_type = factor(count_type, levels=c("unstranded", "sense", "antisense")))
}

#' Create a read count matrix from selected way of counting
read_matrix = function(reads, sel) {
    counts = reads %>% filter(!grepl("^N_", feature)) %>%
        select(sample, feature, {{ sel }}) %>%
        tidyr::pivot_wider(names_from=sample, values_from={{ sel }})

    cmat = data.matrix(counts[-1])
    rownames(cmat) = counts$feature
    cmat[rowSums(cmat) > 0,]
}

#' ggplot2 read count barplot for all samples (mapped, unmapped, ambiguous etc.)
read_barplot = function(stats) {
    ggplot(stats, aes(x=forcats::fct_reorder(sample, -reads), y=reads, fill=feature)) +
        geom_col() +
        facet_wrap(~ count_type, ncol=1, scales="free_y") +
        scale_fill_brewer(palette="Set1") +
        theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'seq_counts/FF230302.rds'),
        opt('p', 'plotfile', 'pdf', 'seq_counts/FF230302.pdf'),
        arg('infiles', '.ReadsPerGene.out.tab', arity='*',
            list.files("seq_aligned/FF230302", "\\.ReadsPerGene\\.out\\.tab", full.names=TRUE))
    )

    reads = read_files(args$infiles)
    stats = read_stats(reads)

    sel = stats %>%
        group_by(count_type) %>%
            summarize(part_mapped = sum(reads[feature == "N_mapped"]) / sum(reads)) %>%
        slice_max(part_mapped) %>%
        pull(count_type)
    levels(stats$count_type)[levels(stats$count_type) == sel] = sprintf("%s (selected)", sel)

    pdf(args$plotfile, 10, 9)
    print(read_barplot(stats))
    dev.off()

    cmat = read_matrix(reads, sel)
    saveRDS(cmat, file=args$outfile)
})
