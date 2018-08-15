library(dplyr)
io = import('io')
seq = import('seq')
sys = import('sys')

#' Return normalized read counts for a count matrix
normalize_reads = function(mat) {
    idx = data.frame(id=colnames(mat))
    DESeq2::DESeqDataSetFromMatrix(mat, colData=idx, ~1) %>%
        DESeq2::estimateSizeFactors() %>%
        DESeq2::counts(normalized=TRUE)
}

extract_segment = function(smp, chr, ratio, genes) {
    `%>%` = magrittr::`%>%`
    message(smp, "", chr)
    use_genes = genes$seqnames == chr
    mat = ratio[use_genes, smp, drop=FALSE]
    ediv = ecp::e.divisive(mat)

    res = cbind(genes[use_genes,], mat=mat, cluster=ediv$cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(start = min(start),
                         end = max(end),
                         ratio = median(mat)) %>%
        dplyr::select(-cluster)

#    density_modal = function(x, bw=1) {
#        x = log2(x[x>0.5])
#        den = density(x, kernel="gaussian", bw=bw)
#        2^den$x[den$y==max(den$y)]
#    }
#    ediv = ecp::e.divisive(as.matrix(data$expr))
#    data$clust = ediv$cluster
#    data %>%
#        group_by(clust) %>%
#        summarize(end = max(start),
#                  start = min(start),
#                  expr = density_modal(expr)) %>%
#        select(-clust)
}

sys$run({
    args = sys$cmd$parse(
        opt('r', 'ref', 'eT expression RData', '../data/rnaseq/Mad2+p53_batch2.RData'),
        opt('e', 'expr', 'expression RData', '../data/rnaseq/assemble.RData'),
        opt('o', 'outfile', 'results RData', 'eT_ploidy.RData'))

    genes = seq$coords$gene("ensembl_gene_id", dset="mmusculus_gene_ensembl", granges=TRUE) %>%
        plyranges::select(ensembl_gene_id) %>%
        as.data.frame() %>%
        filter(seqnames %in% c(1:19, 'X'))

    smps = c("eT_p0", "eT_p2")
    ref = io$load(args$ref)$counts[,smps]
    panel = io$load(args$expr)$counts
    eset = normalize_reads(narray::stack(ref, panel, along=2))
    ref = eset[,colnames(ref)]
    panel = eset[,colnames(panel)]

    keep_fun = function(x) mean(x > 50 & x < 1e4) >= 0.75
    keep_ref = narray::map(ref, along=2, keep_fun)
    keep_panel = narray::map(panel, along=2, keep_fun)
    keep = keep_ref[keep_ref & keep_panel]

    narray::intersect(ref, panel, keep, genes$ensembl_gene_id, along=1)
    ratio = panel / rowMeans(ref)

    segments = expand.grid(sample=colnames(ratio), seqnames=c(1:19,'X'),
                           stringsAsFactors=FALSE) %>%
        mutate(result = clustermq::Q(extract_segment, smp=sample, chr=seqnames,
            const=list(genes=genes, ratio=ratio), n_jobs=15, memory=1024)) %>%
        tidyr::unnest()

    save(segments, genes, file="eT_ploidy.RData")
})
