library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')

# are genes per chr cor w/ number of reads?
chr_df_plot = function(min_reads=10) {
    ctl = readr::read_tsv("../../aneuploidy/data/rnaseq/T-ALL_read_count.txt")
    ctl_reads = data.matrix(ctl[,c("eT_p0","eT_p2")])
    rownames(ctl_reads) = ctl$ensembl_gene_id
    ctl_reads = ctl_reads[rowMeans(ctl_reads) > min_reads,]
    chrs = idmap$gene(rownames(ctl_reads), to="chromosome_name", dset="mmusculus_gene_ensembl")
    chr_genes = data.frame(chr=unname(chrs), genes=names(chrs))
    chr_reads = data.frame(chr=chrs,reads=ctl_reads[,2])
    chrdf = dplyr::inner_join(chr_genes, chr_reads, by="chr") %>%
        group_by(chr) %>%
        summarize(n_genes=n(), n_reads = sum(reads)) %>%
        filter(chr %in% c(1:19, 'X'))

    # calculate correlation stats
    mod = lm(n_reads ~ n_genes, data=chrdf)

    p = ggplot(chrdf, aes(x=n_genes, y=n_reads, label=chr)) +
        geom_text() +
        geom_smooth(method="lm") +
        labs(x = "Number of genes w/ >= 10 reads per sample",
             y = "Number of aligned reads to that chromosome",
             title = sprintf("Correlation of gene number vs expression (p=%.2g, r^2=%.2f)",
                 mod %>% broom::tidy() %>% filter(term == "n_genes") %>% pull(p.value),
                 mod %>% broom::glance() %>% pull(r.squared)))
}

pdf("chr_ngenes_nreads_cor.pdf")
print(chr_df_plot(-1) + labs(subtitle="all genes"))
print(chr_df_plot(10) + labs(subtitle="genes w/ >= 10 reads avg"))
dev.off()

# median, quartiles?
