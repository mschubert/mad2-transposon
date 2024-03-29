library(dplyr)
library(plyranges)
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'transposon insert rds', 'analysis_set.rds'),
    opt('m', 'meta', 'rds', '../data/meta/meta.rds'),
    opt('u', 'upstream', 'bp to include before gene', '10000'),
    opt('d', 'downstream', 'bp to include after gene', '0'),
    opt('o', 'outfile', 'insertion statistics rds', 'poisson.rds'),
    opt('p', 'plotfile', 'volcano pdf', 'poisson.pdf')
)

meta = readRDS(args$meta) %>% filter(analysis_set & !is.na(type))
ins = readRDS(args$infile) %>%
    filter(sample %in% meta$sample,
           chr %in% c(1:19, 'X'),
           !is_local) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,
        start.field="position", end.field="position")

genome = seq$genome("GRCm38", chrs=c(1:19, 'X'))
genes = seq$coords$gene(dset="mmusculus_gene_ensembl", assembly="GRCm38", granges=TRUE) %>%
    filter(gene_biotype == "protein_coding",
           ! grepl("^Gm[0-9]{3,5}", external_gene_name),
           seqnames(.) %in% seqnames(ins)) %>%
    select(external_gene_name) %>%
    group_by(seqnames, strand, external_gene_name) %>% # 100 dups w/ diff pos
    summarize(start = min(start), end = max(end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    anchor_3p() %>% stretch(as.integer(args$upstream)) %>%
    anchor_5p() %>% stretch(as.integer(args$downstream)) %>%
    mutate(TTAAs = seq$count_pattern("TTAA", genome, ., rc=TRUE)) %>%
    filter(TTAAs > 0) # if no upstream, Mir5136 has 0 TTAAs, 1 insert
dists = join_nearest(ins, genes, distance=TRUE)

ins_sites_genome = seq$count_pattern("TTAA", genome, rc=TRUE)
n_smp = length(unique(ins$sample))
sample_rates = as.data.frame(ins) %>%
    group_by(sample) %>%
    summarize(n_ins = dplyr::n(),
              n_reads = sum(reads)) %>%
    mutate(rate = n_ins / ins_sites_genome)
ins_rate_genome = mean(sample_rates$rate)

samples = seq$intersect(ins, genes) %>%
    select(sample, external_gene_name, reads) %>%
    group_by(sample, external_gene_name) %>%
    summarize(n_ins = dplyr::n(),
              reads = sum(reads)) %>%
    ungroup()

ptest = function(n_ins, len, rate) broom::tidy(poisson.test(n_ins, len*n_smp, rate))
result = as.data.frame(GenomicRanges::mcols(genes)) %>%
    inner_join(samples %>% select(sample, external_gene_name)) %>%
    group_by(external_gene_name, TTAAs) %>%
        summarize(n_smp = n_distinct(sample)) %>%
    ungroup() %>%
    mutate(result = purrr::map2(n_smp, TTAAs, ptest, rate=ins_rate_genome)) %>%
    tidyr::unnest() %>%
    select(-(parameter:alternative)) %>%
    mutate(estimate = log2(estimate / ins_rate_genome),
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

p = result %>%
    mutate(size = n_smp,
           label = external_gene_name) %>%
    plt$volcano(label_top=30, repel=TRUE, ceil=1e-30) +
        xlab("log2 fold change poisson rate")

pdf(args$plotfile)
print(p)
dev.off()

saveRDS(list(samples=samples, genes=genes, sample_rates=sample_rates,
             result=result, dists=dists), file=args$outfile)
