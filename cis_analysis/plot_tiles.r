# test for insertion differences with aneuploidy
# plot the result as a volcano (dna, rna separate)/heatmap (both dna+rna)
library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')
sys = import('sys')

plot_dna_tiles = function(cis_samples, cis_result, aneup_assocs, n_dna) {
    genes = cis_result %>%
        filter(adj.p < 1e-5,
               n_smp >= as.integer(n_dna)) %>%
        inner_join(aneup_assocs %>% select(external_gene_name, ast=statistic)) %>%
        arrange(ast) %>%
        pull(external_gene_name)

    dna_tiles = cis_samples %>%
        transmute(sample = sample,
                  gene_name = factor(external_gene_name, levels=genes),
                  ins = 1) %>%
        na.omit() %>%
        tidyr::complete(sample, gene_name, fill=list(ins=0))

    plot_tiles(dna_tiles) +
        theme(legend.position = "left") +
        ggtitle(paste("CIS min", n_dna, "samples, 20 reads"))
}

plot_rna_tiles = function(ins_rna, ctg, valid_sample, n_rna) {
    rna_tiles = ins_rna %>%
        select(sample, gene_name) %>%
        mutate(ins = 1) %>%
        tidyr::complete(sample, gene_name, fill=list(ins=0)) %>%
        filter(sample %in% valid_sample,
               gene_name %in% ctg$gene_name) %>%
        na.omit() %>%
        tidyr::complete(sample, gene_name, fill=list(ins=NA))

    plot_tiles(rna_tiles) +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank()) +
        guides(fill=FALSE) +
        ggtitle(paste("CTGs min", n_rna, "samples"))
}

plot_tiles = function(data) {
    ggplot(data, aes(x=gene_name, y=sample)) +
        geom_tile(aes(fill=ins), color="white") +
        coord_fixed() +
        viridis::scale_fill_viridis(option="magma", direction=-1, na.value="#f5f5f5") +
        theme(axis.text.x = element_text(size=10, angle=65, hjust=1),
              axis.text.y = element_text(size=10),
              axis.title.x = element_text(size=12))
}

plot_aneup = function(aneup) {
    ggplot(aneup, aes(x=aneuploidy, y=sample)) +
        geom_segment(aes(xend=aneuploidy, yend=sample), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.x = element_text(size=10),
              axis.title.x = element_text(size=12),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aneup', 'sample-level aneup scores', '../ploidy_compare/analysis_set.RData'),
        opt('i', 'ins_dna', 'all DNA insertions table', 'poisson.RData'),
        opt('j', 'ins_rna', 'all RNA insertions table', '../data/rnaseq_imfusion/insertions.txt'),
        opt('d', 'assocs_dna', 'CIS in DNA', 'aneup_assocs.RData'),
        opt('r', 'assocs_rna', 'CTG in RNA', '../data/rnaseq_imfusion/merged_ctgs.txt'),
        opt('e', 'exons', 'exon expression table', '../data/rnaseq_imfusion/exon_counts.txt'),
        opt('n', 'n_dna', 'min number samples with ins', '5'),
        opt('m', 'n_rna', 'min number samples with ins', '2'),
        opt('p', 'plotfile', 'pdf to plot to', 'plot_tiles.pdf'))

    aneup_assocs = io$load(args$assocs_dna)$aneuploidy
    aneup = io$load(args$aneup) %>%
        arrange(aneuploidy) %>%
        mutate(sample = factor(sample, levels=sample))

    cis = io$load(args$ins_dna)
    cis_result = cis$result
    cis_samples = cis$samples %>%
        mutate(sample = factor(sample, levels=levels(aneup$sample)))

    expr = io$read_table(args$exons, header=TRUE)
    ins_rna = io$read_table(args$ins_rna, header=TRUE) %>%
        mutate(sample = factor(sample, levels=levels(aneup$sample)))
    ctg = io$read_table(args$assocs_rna, header=TRUE) %>%
        filter(n_samples >= as.integer(args$n_rna))

    left = plot_dna_tiles(cis_samples, cis_result, aneup_assocs, args$n_dna)
    mid = plot_rna_tiles(ins_rna, ctg, colnames(expr), args$n_rna)
    right = plot_aneup(aneup)

    pdf(args$plotfile, 16, 14)
    print(left + mid + right + plot_layout(nrow=1, widths=c(5,2,2)))
    dev.off()
})
