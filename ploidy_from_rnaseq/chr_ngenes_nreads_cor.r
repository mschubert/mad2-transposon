library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')

#' See if gene number correlates w/ reads across chromosomes
#'
#' Reason for this is that for ploidy estimates from RNA-seq, it would be good
#' if we could fit one parameter per sample, instead of one parameter per sample
#' per chromosome. This looks like the following:
#'
#'    n_reads
#' --------------  =  const per sample
#' somy * n_genes
#'
#'       n_reads^chr
#' ----------------------  =  k^genome * n_reads^genome
#' somy^chr * n_genes^chr
#'
#' or we should rather go by chromosome:
#'
#'                   n_reads^chr
#' somy = k^chr *  -------------- = k^chr * frac_reads
#'                 n_reads^genome
#'
#' This data indicates that we should adjust per chromosome, as the correlation
#' is not stable enough across samples (but better then chromosome length).
#'
#' So use single-cell WGS to fit k^chr for each chromosome, and then use those
#' to infer ploidy for the samples where we don't have it.
invisible(NULL)

#' Summarize number of reads and genes per chromosome
#'
#' @param reads      Matrix of genes (rows) x samples (columns)
#' @param min_reads  How many reads need to be there for a gene to be considered
#' @return           Data.frame with fields: n_genes, n_reads, frac_reads
chr_read_genes = function(reads, min_reads=10, chromosomes=c(1:19, 'X'),
                       dset="mmusculus_gene_ensembl", facet=TRUE) {
    chrs = idmap$gene(rownames(reads), to="chromosome_name", dset=dset)
    rownames(reads) = unname(chrs)
    names(dimnames(reads)) = c("chr", "sample")
    read_df = reshape2::melt(reads, value.name="n_reads") %>%
        filter(n_reads >= min_reads,
               chr %in% chromosomes) %>%
        group_by(sample) %>%
        mutate(frac_reads = n_reads / sum(n_reads)) %>%
        ungroup() %>%
        group_by(chr, sample) %>%
        summarize(n_genes = n(),
                  n_reads = sum(n_reads),
                  frac_reads = sum(frac_reads))
}

#' Plot the correlation between number of reads and number of genes per chromosome
#'
#' @param read_df  Data.frame from chr_read_genes
#' @return         ggplot2 object
plot_all = function(read_df) {
    mod = lm(n_reads ~ n_genes, data=read_df)

    p = ggplot(read_df, aes(x=n_genes, y=frac_reads, label=chr, color=sample)) +
        geom_text() +
        geom_smooth(method="lm", se=FALSE) +
        labs(x = "Number of expressed genes",
             y = "Fraction of aligned reads to that chromosome",
             title = sprintf("Correlation of gene number vs expression (p=%.2g, r^2=%.2f)",
                 mod %>% broom::tidy() %>% filter(term == "n_genes") %>% pull(p.value),
                 mod %>% broom::glance() %>% pull(r.squared)))
}

#' Plot the correlation between number of reads and number of genes per chromosome
#'
#' @param read_df  Data.frame from chr_read_genes
#' @return         Facetted ggplot2 object
plot_facets = function(read_df) {
    p = ggplot(read_df, aes(x=n_genes, y=frac_reads, label=chr)) +
        geom_text() +
        geom_smooth(method="lm") +
        facet_wrap(~ sample)
}

if (is.null(module_name())) {
    tumors = readr::read_tsv("../../aneuploidy/data/rnaseq/T-ALL_read_count.txt")
    tumor_reads = data.matrix(tumors[,c("eT_p0","eT_p2")])
    rownames(tumor_reads) = tumors$ensembl_gene_id

    normals = readr::read_tsv("../../aneuploidy/data/rnaseq/seq_run2.txt")
    normal_reads = data.matrix(normals[,c("lib11", "lib12")])
    colnames(normal_reads) = c("T372 (wt)", "T373 (wt)")
    rownames(normal_reads) = normals$Gene

    reads = narray::stack(list(tumor_reads, normal_reads), along=2, fill=0)

    dfs = list(`0` = chr_read_genes(reads, 0),
               `10` = chr_read_genes(reads, 10),
               `500` = chr_read_genes(reads, 500))

    pdf("chr_ngenes_nreads_cor.pdf")
    for (d in seq_along(dfs)) {
        p1 = plot_all(dfs[[d]])
        p2 = plot_facets(dfs[[d]])
        subt = labs(subtitle = sprintf("genes w/ >= %s reads", names(dfs)[d]))
        print(p1 + subt)
        print(p2 + subt)
    }
    dev.off()
}
