library(dplyr)
seq = import('seq')
sys = import('sys')

#' Return normalized read counts for a count matrix
#'
#' @param mat  gene expression matrix [genes x samples]
#' @return  normalized read counts
normalize_reads = function(mat) {
    idx = data.frame(id=colnames(mat))
    DESeq2::DESeqDataSetFromMatrix(mat, colData=idx, ~1) %>%
        DESeq2::estimateSizeFactors() %>%
        DESeq2::counts(normalized=TRUE)
}

#' Center a vector by its density maximum
#'
#' @param x  numeric vector
#' @param bw  bin width for the density estimator
#' @return  vector centered by its density maximum
center_segment_density = function(x, w=NULL, bw=bw.nrd0(x)) {
    den = density(x, kernel="gaussian", bw=bw, weights=w)
    x - den$x[den$y==max(den$y)]
}

#' Calculate continuous copy segments along chromosome
#'
#' @param smp  character string for sample identifier
#' @param chr  character string for chromosome
#' @param ratio  numeric matrix [genes x samples]
#' @param genes  gene info for rows in ratio
#' @param bw  kernel width for estimating density (in ploidies up/down)
#' @return  data.frame with estimated ploidy segments
extract_segment = function(smp, chr, ratio, genes, bw=NULL) {
    center_of_density = function(x, bw) {
        if (is.null(bw))
            bw = bw.nrd0(x)
#        x = log2(x[x>0.5])
        den = density(x, kernel="gaussian", bw=bw)
        den$x[den$y==max(den$y)]
#        2^den$x[den$y==max(den$y)]
    }

    `%>%` = magrittr::`%>%`
    message(smp, "", chr)
    chr_genes = genes %>%
        dplyr::filter(seqnames == chr) %>%
        dplyr::arrange(start)
    mat = ratio[chr_genes$ensembl_gene_id, smp]
    ediv = ecp::e.divisive(as.matrix(mat), min.size=50)

    res = cbind(chr_genes, cmat=mat, cluster=ediv$cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(start = min(start),
                         end = max(end),
                         width = abs(end - start),
                         ploidy = 2 * center_of_density(cmat, bw=bw)) %>%
        dplyr::select(-cluster)
}

sys$run({
    args = sys$cmd$parse(
        opt('r', 'ref', 'eT expression rds', '../data/rnaseq/Mad2+p53_batch2.rds'),
        opt('e', 'expr', 'expression rds', '../data/rnaseq/assemble.rds'),
        opt('o', 'outfile', 'results rds', 'eT_ploidy.rds'))

    genes = seq$coords$gene("ensembl_gene_id", dset="mmusculus_gene_ensembl", granges=TRUE) %>%
        plyranges::select(ensembl_gene_id) %>%
        as.data.frame() %>%
        filter(seqnames %in% c(1:19, 'X'))

    smps = c("eT_p0", "eT_p2")
    ref = readRDS(args$ref)$counts[,smps]
    panel = readRDS(args$expr)$counts
    eset = normalize_reads(narray::stack(ref, panel, along=2))
    ref = eset[,colnames(ref)]
    panel = eset[,colnames(panel)]

    keep_fun = function(x) mean(x > 50 & x < 1e4) >= 0.75
    keep_ref = narray::map(ref, along=2, keep_fun)
    keep_panel = narray::map(panel, along=2, keep_fun)
    keep = keep_ref[keep_ref & keep_panel]

    narray::intersect(ref, panel, keep, genes$ensembl_gene_id, along=1)
    ratio = panel / rowMeans(ref)

    segments = expand.grid(sample=colnames(ratio), seqnames=c(1:19,'X')) %>%
        mutate(sample = as.character(sample),
               result = clustermq::Q(extract_segment, smp=sample, chr=seqnames,
                   const=list(genes=genes, ratio=ratio),
                   n_jobs=15, memory=1024)) %>%
        tidyr::unnest() %>%
        group_by(sample) %>%
        mutate(ploidy = 2 + center_segment_density(ploidy, w=width, bw=0.25)) %>%
        ungroup()

    ratio = cbind(genes, ratio) %>%
        tidyr::gather("sample", "ratio", -(seqnames:ensembl_gene_id))

    saveRDS(list(segments=segments, ratio=ratio), file="eT_ploidy.rds")
})
