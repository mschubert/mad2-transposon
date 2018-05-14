# test for insertion differences with aneuploidy
# plot the result as a volcano (dna, rna separate)/heatmap (both dna+rna)
library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')
sys = import('sys')

#' Plot read counts for all samples/exons with anntotation insertion sites
plot_exon_expr = function() {
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aneuploidy', 'sample-level aneup scores', '../tis2_assoc-tryout/dset.RData'),
        opt('i', 'ins_dna', 'all DNA insertions table', 'cis_sanger.tsv'),
        opt('j', 'ins_rna', 'all RNA insertions table', '../data/rnaseq_imfusion/insertions.txt'),
        opt('d', 'assocs_dna', 'CIS in DNA', 'cis_sanger_results.tsv'),
        opt('r', 'assocs_rna', 'CTG in RNA', '../data/rnaseq_imfusion/merged_ctgs.txt'),
        opt('e', 'exons', 'exon expression table', '../data/rnaseq_imfusion/exon_counts.txt'),
        opt('p', 'plotfile', 'pdf to plot to', 'cis_tiles.pdf'))

    expr = io$read_table(args$exons, header=TRUE)

    aneup = io$load(args$aneuploidy) %>%
        select(sample, aneup) %>%
        distinct() %>%
        na.omit() %>%
        arrange(aneup) %>%
        mutate(sample = factor(sample, levels=sample))

#    ins_dna = io$read_table(args$ins_dna, header=TRUE)
    ins_dna = io$load(args$aneuploidy) %>%
        filter(reads >= 20) %>% select(-aneup) # has the right sample identifiers (fix?)
    cis = io$read_table(args$assocs_dna, header=TRUE)
    ks = io$load("../tis2_assoc-tryout/ks.RData")$hits %>%
        filter(ensembl_gene_id %in% cis$gene_id) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(statistic)
    dna_tiles = ins_dna %>%
        filter(gene_name %in% ks$gene_name) %>%
        transmute(sample = factor(sample, levels=levels(aneup$sample)),
                  gene_name = factor(gene_name, levels=ks$gene_name),
                  ins = 1) %>%
        na.omit() %>%
        tidyr::complete(sample, gene_name, fill=list(ins=0))
    left = ggplot(dna_tiles, aes(x=gene_name, y=sample)) +
        geom_tile(aes(fill=ins), color="white") +
        coord_fixed() +
        viridis::scale_fill_viridis(option="magma", direction=-1) +
        theme(axis.text.x = element_text(size=10, angle=65, hjust=1),
              axis.text.y = element_text(size=10),
              axis.title.x = element_text(size=12),
              legend.position = "left") +
        ggtitle("CIS min 5 samples, 20 reads")

    ins_rna = io$read_table(args$ins_rna, header=TRUE)
    ctg = io$read_table(args$assocs_rna, header=TRUE) %>%
        filter(n_samples > 1)
    rna_tiles = ins_rna %>%
        transmute(sample = factor(sample, levels=grep("^[0-9]+[st]$", colnames(expr), value=TRUE)), # all RNA
                  gene_name = gene_name,
                  ins = 1) %>%
        tidyr::complete(sample, gene_name, fill=list(ins=0)) %>%
        filter(gene_name %in% ctg$gene_name) %>%
        mutate(sample = factor(sample, levels=levels(aneup$sample))) %>%
        na.omit() %>%
        tidyr::complete(sample, gene_name, fill=list(ins=NA))
    mid = ggplot(rna_tiles, aes(x=gene_name, y=sample)) +
        geom_tile(aes(fill=ins), color="white") +
        coord_fixed() +
        viridis::scale_fill_viridis(option="magma", direction=-1, na.value="transparent") +
        theme(axis.text.x = element_text(size=10, angle=65, hjust=1),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size=12)) +
        guides(fill=FALSE) +
        ggtitle("CTGs min 2 samples")

    right = ggplot(aneup, aes(x=aneup, y=sample)) +
        geom_segment(aes(xend=aneup, yend=sample), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.x = element_text(size=10),
              axis.title.x = element_text(size=12),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())

    pdf(args$plotfile, 10, 8)
    print(left + mid + right + plot_layout(nrow=1, widths=c(3,1,3)))
    dev.off()

})
