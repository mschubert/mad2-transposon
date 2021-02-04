library(dplyr)
library(plyranges)
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'transposon insert RData', 'analysis_set.RData'),
    opt('f', 'feature', 'enhancer|promoter|TF_binding_site etc.', 'enhancer'),
    opt('u', 'upstream', 'bp to include before gene', '5000'),
    opt('d', 'downstream', 'bp to include after gene', '5000'),
    opt('o', 'outfile', 'insertion statistics RData', 'poisson_epi.RData'),
    opt('p', 'plotfile', 'volcano pdf', 'poisson_epi.pdf'))

ins = io$load(args$infile) %>%
    filter(chr %in% c(1:19, 'X'), !is_local) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,
        start.field="position", end.field="position")

genome = seq$genome("GRCm38", chrs=c(1:19, 'X'))

genes = seq$coords$gene(dset="mmusculus_gene_ensembl", granges=TRUE)
# time out
## biomart query ensembl mouse regulation
#ensembl = biomaRt::useEnsembl("ENSEMBL_MART_FUNCGEN", dataset="mmusculus_regulatory_feature")
#attrs = c("chromosome_name", "chromosome_start", "chromosome_end", "feature_type_name",
#          "regulatory_stable_id", "epigenome_name", "activity")
#res = biomaRt::getBM(attributes=attrs, filters=c("epigenome_name", "activity"),
#                     values=list(c("spleen adult","thymus adult"), "ACTIVE"), mart=ensembl)
res = readr::read_tsv("poisson_epi.tsv") %>%
    filter(`SO term name` == args$feature) %>%
    dplyr::rename(chr = `Chromosome/scaffold name`,
                  start = `Start (bp)`,
                  end = `End (bp)`,
                  id = `Regulatory stable ID`) %>%
    select(chr, start, end, id) %>%
    distinct() %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    filter(seqnames(.) %in% seqnames(ins)) %>%
    anchor_center %>% stretch(10000) %>%
    mutate(TTAAs = seq$count_pattern("TTAA", genome, ., rc=TRUE)) %>%
    filter(TTAAs > 0) # if no upstream, Mir5136 has 0 TTAAs, 1 insert
res = res[count_overlaps(res, genes) == 0]

ins_sites_genome = seq$count_pattern("TTAA", genome, rc=TRUE)
n_smp = length(unique(ins$sample))
sample_rates = as.data.frame(ins) %>%
    group_by(sample) %>%
    summarize(n_ins = dplyr::n(),
              n_reads = sum(reads)) %>%
    mutate(rate = n_ins / ins_sites_genome)
ins_rate_genome = mean(sample_rates$rate)

samples = seq$intersect(ins, res) %>%
    select(sample, id, reads) %>%
    group_by(sample, id) %>%
    summarize(n_ins = dplyr::n(),
              reads = sum(reads)) %>%
    ungroup()

ptest = function(n_ins, len, rate) broom::tidy(poisson.test(n_ins, len*n_smp, rate))
result = as.data.frame(GenomicRanges::mcols(res)) %>%
    inner_join(samples %>% select(sample, id)) %>%
    group_by(id, TTAAs) %>%
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
           label = id) %>%
    plt$p_effect() %>%
    plt$volcano(label_top=30, repel=TRUE, ceil=1e-30) +
        xlab("log2 fold change poisson rate")

pdf(args$plotfile)
print(p)
dev.off()

save(samples, res, sample_rates, result, file=args$outfile)
