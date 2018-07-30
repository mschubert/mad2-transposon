library(dplyr)
library(plyranges)
io = import('io')
seq = import('seq')
sys = import('sys')

#TODO: only allow one insertion per sample and region?

args = sys$cmd$parse(
    opt('i', 'infile', 'cis RData', 'analysis_set.RData'),
    opt('u', 'upstream', 'bp to include before gene', '1000'), #TOOD: 10k (?)
    opt('d', 'downstream', 'bp to include after gene', '0'),
    opt('o', 'outfile', 'insertion statistics RData', 'poisson.RData'))

cis = io$load(args$infile) %>%
    filter(chr %in% c(1:19, 'X')) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,
        start.field="position", end.field="position")

genome = seq$genome("GRCm38", chrs=c(1:19, 'X'))
ins_sites_total = seq$count_pattern("TTAA", genome, rc=TRUE)
ins_per_sample = length(cis) / length(unique(cis$sample))
ins_rate = ins_per_sample / ins_sites_total

#stretch(anchor_3p(genes), as.integer(args$upstream)) #TODO: check if this works
#stretch(anchor_5p(genes), as.integer(args$downstream))
genes = seq$coords$gene(dset="mmusculus_gene_ensembl", granges=TRUE) %>%
    filter(seqnames(.) %in% seqnames(cis)) %>%
    select(external_gene_name) %>% # ...stretch here...
    mutate(TTAAs = seq$count_pattern("TTAA", genome, ., rc=TRUE),
           cis = count_overlaps(., cis)) %>%
    arrange(-cis)

do_test = function(n_ins, gene_len, rate)
    poisson.test(n_ins, gene_len * n_smp, rate)$p.value

result = tbl_df(as.data.frame(genes)) %>%
    mutate(p.value = purrr::map2_dbl(cis, TTAAs, do_test, rate=ins_rate),
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

save(result, file=args$outfile)
