library(reshape2)
library(dplyr)
library(ggplot2)
library(ggbio)
library(GenomicRanges)
io = import('io')
seq = import('seq')

ref = io$load("../data/rnaseq/Mad2+p53_batch2.RData")$counts[,c("eT_p0", "eT_p2")]
ctl = edgeR::cpm(ref)

load('../data/rnaseq/assemble.RData') # counts, expr, idx
cpm = edgeR::cpm(counts)
narray::intersect(ctl, ref, cpm, counts, along=1)
keep = narray::map(cbind(ctl, cpm), along=2, function(x) all(x > 0.5 & x < 1000))
ctl = edgeR::cpm(ref[keep,])
cpm = edgeR::cpm(counts[keep,])
scaled = (2 * cpm / narray::crep(rowMeans(ctl), ncol(cpm))) %>%
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
    density_modal = function(x, bw=1) {
        x = log2(x[x>0.5])
        den = density(x, kernel="gaussian", bw=bw)
        2^den$x[den$y==max(den$y)]
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
    mutate(sample = sub("^X", "", sample)) %>%
    tbl_df()

plot_sample = function(smp) {
    cur = segments %>% filter(sample == smp)
    autoplot(genes, aes_string(y=smp), geom="point", shape=1, alpha=0.3) +
        geom_segment(data=cur, aes(x=start, xend=end, y=expr, yend=expr), size=3, color="green") +
        scale_y_continuous(trans="log2", breaks=c(1:6)) +
#        geom_hline(aes(yintercept=median(genes[[smp]])), linetype="dashed", color="grey", size=2) +
        coord_cartesian(ylim=c(0.5,6)) +
        theme(axis.text.x = element_blank()) +
        ggtitle(smp)
}

pdf("eT_ploidy.pdf", 10, 4)
for (smp in unique(segments$sample))
    print(plot_sample(smp))
dev.off()

save(segments, genes, file="eT_ploidy.RData")
