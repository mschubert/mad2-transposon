library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

plot_pca = function(eset) {
    rc = colSums(DESeq2::counts(eset, normalized=FALSE))
    vst = DESeq2::varianceStabilizingTransformation(eset[,rc>0])
    pcadata = DESeq2::plotPCA(vst, intgroup=c("sample_id", "mouse", "cell_type", "treatment"), returnData=TRUE) %>%
        mutate(sample_id = sub(".*_([0-9]+)$", "\\1", sample_id))
    pcavar = round(100 * attr(pcadata, "percentVar"))

    ggplot(pcadata, aes(PC1, PC2)) +
        geom_point(aes(shape=treatment, color=cell_type), size=3, alpha=0.8) +
        ggrepel::geom_text_repel(aes(label=sample_id), max.overlaps=Inf) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance"))
}

file2deseq = function(fname, samples) {
    rtab = readr::read_tsv(fname)
    reads = data.matrix(rtab[-1])
    rownames(reads) = rtab$Gene
    reads = reads[,samples$sample_id]

    eset = DESeq2::DESeqDataSetFromMatrix(reads, samples, ~1) #%>%
#        DESeq2::estimateSizeFactors()
}

plot_overview = function(eset, title="") {
    rc = tibble(sample_id = sub(".*_([0-9]+)$", "\\1", colnames(eset)),
                read_count = colSums(DESeq2::counts(eset, normalized=FALSE)),
                gene_count = colSums(DESeq2::counts(eset, normalized=FALSE) != 0)) %>%
        mutate(sample_id = factor(sample_id, levels=gtools::mixedsort(sample_id)))
    prc = ggplot(rc, aes(x=sample_id, y=read_count)) + geom_col() + ggtitle(title)
    pgc = ggplot(rc, aes(x=sample_id, y=gene_count)) + geom_col()
    ppca = plot_pca(eset)

    (prc / pgc) | ppca
}

sys$run({
    args = sys$cmd$parse(
        opt('s', 'samples', 'tsv', 'samples.tsv'),
        opt('o', 'outfile', 'rds', 'rnaseq.rds'),
        opt('p', 'plotfile', 'pdf', 'rnaseq.pdf')
    )

    tab = readxl::read_excel("Total_gene_counts.xls")
    samples = readr::read_tsv(args$samples)

    fnames = c("SCRB_CH210203_UMIcounts_only_plus.txt.gz",
               "SCRB_CH210203_UMIcounts_only_minus.txt.gz",
               "SCRB_CH210203_UMIcounts_both_strands.txt.gz")

    esets = lapply(fnames, file2deseq, samples=samples)
    plots = mapply(plot_overview, eset=esets, title=fnames, SIMPLIFY=FALSE)

    pdf(args$plotfile, 12, 7)
    print(ggplot() + annotation_custom(gridExtra::tableGrob(samples)))
    for (p in plots)
        print(p)
    dev.off()
})
