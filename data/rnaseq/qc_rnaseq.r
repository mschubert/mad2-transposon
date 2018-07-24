library(dplyr)
library(ggrepel)
library(RColorBrewer)
b = import('base')
io = import('io')
sys = import('sys')
plt = import('plot')
rnaseq = import('process/rna-seq')

#' do PCA/dim reduction plots to see how they cluster
plot_pca = function(expr, idx) {
    if (is.null(idx$tissue))
        idx$tissue = "unknown"
    if (is.null(idx$type))
        idx$type = "unknown"

    pca = prcomp(t(expr), scale=FALSE)
    p1 = ggplot(cbind(idx, pca$x), aes(x=PC1, y=PC2, color=tissue, shape=type)) +
        geom_point(size=5) +
        geom_text_repel(aes(label=sample), color="black") +
        labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
             y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
             title = "PCA plot (linear)")
}

#' same for tsne
plot_tsne = function(expr, idx) {
    if (is.null(idx$tissue))
        idx$tissue = "unknown"
    if (is.null(idx$type))
        idx$type = "unknown"

    tsne = Rtsne::Rtsne(t(expr), perplexity=round(ncol(expr)/3)-1)
    p2 = cbind(idx, x=tsne$Y[,1], y=tsne$Y[,2]) %>%
        ggplot(aes(x=x, y=y, color=tissue, shape=type)) +
        geom_point(size=5) +
        geom_text_repel(aes(label=sample), color="black") +
        labs(x = "tsne 1",
             y = "tsne 2",
             title = "T-SNE plot (non-linear)")
}

#' Plot read count statistics
plot_stats = function(stats) {
    colnames(stats) = tools::file_path_sans_ext(basename(colnames(stats)))
    stats = stats %>%
        tidyr::gather("sample", "reads", -Status) %>%
        filter(reads != 0) %>%
        mutate(sample = factor(sample),
               Status = factor(Status, levels=rev(unique(sort(Status)))))
    p = ggplot(stats, aes(x=reorder(sample, -reads, sum), y=reads, fill=Status)) +
        geom_bar(stat="identity") +
        geom_hline(yintercept=5e6, linetype="dashed") +
        labs(x = "sample",
             title = sprintf("%.2f M reads total (avg %.2f M per sample)",
                             sum(stats$reads)/1e6,
                             sum(stats$reads)/1e6/length(unique(stats$sample)))) +
        theme(axis.text.x = element_text(angle=45, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'read count tsv', 'Mad2+PB_batch3.tsv'),
        opt('s', 'stats', 'STAR summary file', 'Mad2+PB_batch3.tsv.summary'),
        opt('m', 'meta', 'sample metadata', '../meta/180620 Overview table mice transposon screen.tsv'),
        opt('o', 'outfile', 'save results to .RData', 'merge_rnaseq.RData'),
        opt('p', 'plotfile', 'qc plot pdf', 'merge_rnaseq.pdf'))

    stats = io$read_table(args$stats, header=TRUE)

    edf = io$read_table(args$infile, header=TRUE, skip=1)
    counts = data.matrix(edf[,-(1:6)])
    colnames(counts) = tools::file_path_sans_ext(basename(colnames(counts)))
    rownames(counts) = edf$Geneid

    if (grepl("PB", args$infile)) {
        idx = io$read_table(args$meta, header=TRUE) %>%
            transmute(hist_nr = `Hist nr.`,
                      type = `Diagnosis (Based on pathology & FACS) (t=thymus, s=spleen, l=liver, ln=lymph node, k=kidney, i=intestine)`,
                      tissue = `RNA sequencing`) %>%
            mutate(tissue = tolower(gsub("[^ST]", "", tissue)),
                   tissue = sapply(tissue, strsplit, split=NULL)) %>%
            tidyr::unnest() %>%
            mutate(type = sapply(type, function(t) b$grep("^([[^ ^/]]+ [[^ ^/]]+)", t)),
                   sample = paste0(hist_nr, tissue))
    } else
        idx = data.frame(sample=colnames(counts))

    narray::intersect(idx$sample, counts, along=2)
    expr = rnaseq$vst(counts)
    dd = as.matrix(dist(t(expr)))

    pdf(args$plotfile, width=10, height=8)
    print(plot_stats(stats))
    pheatmap::pheatmap(dd, col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
    print(plot_pca(expr, idx))
    print(plot_tsne(expr, idx))
    dev.off()

    save(expr, counts, idx, file=args$outfile)
})
