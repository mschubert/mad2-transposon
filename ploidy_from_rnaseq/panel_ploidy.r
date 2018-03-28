library(reshape2)
library(dplyr)
library(ggplot2)
library(ggbio)
library(GenomicRanges)
io = import('io')
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
    as.data.frame()
colnames(scaled) = paste0("x", colnames(scaled))
scaled=scaled%>%
    mutate(ensembl_gene_id = rownames(.))

coords = seq$coords$gene("ensembl_gene_id", dset="mmusculus_gene_ensembl", granges=TRUE) %>%
    plyranges::select(ensembl_gene_id) %>%
    as.data.frame() %>%
    dplyr::inner_join(scaled) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
genome = seq$genome("GRCm38")
seqinfo(coords) = seqinfo(genome)

extract_segment = function(data) {
    ediv = ecp::e.divisive(as.matrix(data$expr))
    data$clust = ediv$cluster
    data %>%
        group_by(clust) %>%
        summarize(start2 = min(start),
                  end = max(start),
                  expr = median(expr)) %>%
        dplyr::rename(start=start2) %>%
        select(-clust)
}
segs = as.data.frame(coords) %>%
    filter(seqnames %in% c(1:19,'X')) %>%
    arrange(seqnames, start) %>%
    tidyr::gather("sample", "expr", -(seqnames:ensembl_gene_id)) %>%
    group_by(seqnames, sample) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, extract_segment)) %>%
    select(-data) %>%
    tidyr::unnest()

scale_modal_diploid = function(x) {
    dens = density(x, kernel="gaussian", bw=(max(x)-min(x))/5)
    2 * x / dens$x[dens$y==max(dens$y)]
}
aneup = segs %>%
#    group_by(sample) %>%
#    mutate(expr = scale_modal_diploid(expr)) %>%
#    ungroup() %>%
    mutate(width = abs(end-start)) %>%
    seq$aneuploidy(sample="sample", ploidy="expr")

plot_sample = function(smp) {
    cur = segs %>% filter(sample == smp)
    autoplot(coords, aes_string(y=smp), geom="point", shape=1, alpha=0.3) +
        geom_segment(data=cur, aes(x=start, xend=end, y=expr, yend=expr), size=3, color="green") +
        scale_y_continuous(trans="log2", breaks=c(1:6)) +
#        geom_hline(aes(yintercept=median(coords[[smp]])), linetype="dashed", color="grey", size=2) +
        coord_cartesian(ylim=c(0.5,6)) +
        theme(axis.text.x = element_blank()) +
        ggtitle(smp)
}

pdf("panel_ploidy.pdf", 10, 4)
for (smp in unique(segs$sample))
    print(plot_sample(smp))
dev.off()

save(aneup, segs, coords)
