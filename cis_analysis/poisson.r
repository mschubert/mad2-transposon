library(dplyr)
library(plyranges)
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'transposon insert RData', 'analysis_set.RData'),
    opt('u', 'upstream', 'bp to include before gene', '10000'),
    opt('d', 'downstream', 'bp to include after gene', '0'),
    opt('o', 'outfile', 'insertion statistics RData', 'poisson.RData'))

ins = io$load(args$infile) %>%
    filter(chr %in% c(1:19, 'X')) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,
        start.field="position", end.field="position")

genome = seq$genome("GRCm38", chrs=c(1:19, 'X'))
ins_sites_total = seq$count_pattern("TTAA", genome, rc=TRUE)
n_smp = length(unique(ins$sample))
ins_per_sample = length(ins) / n_smp
ins_rate = ins_per_sample / ins_sites_total

genes = seq$coords$gene(dset="mmusculus_gene_ensembl", granges=TRUE) %>%
    filter(seqnames(.) %in% seqnames(ins)) %>%
    select(external_gene_name) %>%
    anchor_3p() %>% stretch(as.integer(args$upstream)) %>%
    anchor_5p() %>% stretch(as.integer(args$downstream)) %>%
    mutate(TTAAs = seq$count_pattern("TTAA", genome, ., rc=TRUE))

samples = seq$intersect(ins, genes) %>%
    select(sample, external_gene_name, reads) %>%
    group_by(sample, external_gene_name) %>%
    summarize(n_ins = dplyr::n(),
              reads = sum(reads)) %>%
    ungroup()

ptest = function(n_ins, len, rate) poisson.test(n_ins, len * n_smp, rate)$p.value
result = as.data.frame(GenomicRanges::mcols(genes)) %>%
    inner_join(samples %>% select(sample, external_gene_name)) %>%
    group_by(external_gene_name, TTAAs) %>%
    summarize(n_smp = n_distinct(sample)) %>%
    mutate(p.value = purrr::map2_dbl(n_smp, TTAAs, ptest, rate=ins_rate),
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

save(samples, genes, result, file=args$outfile)
