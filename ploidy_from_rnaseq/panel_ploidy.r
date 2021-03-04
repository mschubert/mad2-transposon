library(reshape2)
library(dplyr)
library(ggplot2)
library(ggbio)
library(GenomicRanges)
seq = import('seq')

load('../data/rnaseq/assemble.RData') # counts, expr, idx
cpm = edgeR::cpm(counts)
cpm = edgeR::cpm(counts[narray::map(cpm, along=2, function(x) all(x > 1 & x < 500)),])
scale_cpm = function(x) {
#    dens = density(x, kernel="gaussian", bw=(max(x)-min(x))/5)
#    log2(x / dens$x[dens$y==max(dens$y)])
    2*x/median(x)
}
scaled = narray::map(cpm, along=2, scale_cpm, subsets=as.character(grepl("[sS]", colnames(cpm)))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id")

genes = seq$coords$gene("ensembl_gene_id", dset="mmusculus_gene_ensembl", granges=TRUE) %>%
    plyranges::select(ensembl_gene_id) %>%
    as.data.frame() %>%
    dplyr::inner_join(scaled) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
genome = seq$genome("GRCm38")
seqinfo(genes) = seqinfo(genome)

extract_segment = function(data) {
    density_modal = function(x) {
        den = density(x, kernel="gaussian", bw=5)
        den$x[den$y==max(den$y)]
    }
    ediv = ecp::e.divisive(as.matrix(data$expr))
    data$clust = ediv$cluster
    data %>%
        group_by(clust) %>%
        summarize(end = max(start),
                  start = min(start),
                  expr = density_modal(expr)) %>%
        select(-clust)
}
segments = as.data.frame(genes) %>%
    filter(seqnames %in% c(1:19,'X')) %>%
    arrange(seqnames, start) %>%
    tidyr::gather("sample", "expr", -(seqnames:ensembl_gene_id)) %>%
    group_by(seqnames, sample) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, extract_segment)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    as.data.frame() %>%
    tbl_df()

saveRDS(list(segments=segments, genes=genes), file="panel_ploidy.rds")
