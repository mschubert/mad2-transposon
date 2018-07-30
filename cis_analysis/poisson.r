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

# collapse close insertion sites per sample?
cis = io$load(args$infile) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,
        start.field="position", end.field="position")

#TODO:
# function: for each chromosome, count TTAA and AATT; also count per gene
# use on all chrs
genome = seq$genome("GRCm38")
# subset default chromosomes
# construct GRanges with all TTAA sites
# count sites in genes and use as gene length
ttaa = Biostrings::matchPattern("TTAA", genome$`1`)
aatt = Biostrings::matchPattern("AATT", genome$`1`)
# count each region only once per sample?

#stretch(anchor_3p(genes), as.integer(args$upstream)) #TODO: check if this works
#stretch(anchor_5p(genes), as.integer(args$downstream))
genes = seq$coords$gene(dset="mmusculus_gene_ensembl", granges=TRUE) %>%
    filter(seqnames(.) %in% seqnames(cis)) %>%
    select(external_gene_name) %>% # ...stretch here...
    mutate(cis = count_overlaps(., cis)) %>%
    arrange(-cis)

glen = seq$lengths(seq$genome("GRCm38"))
glen = glen[names(glen) %in% seqnames(cis)]
n_smp = length(unique(cis$sample))
rate = sum(genes$cis) / (sum(glen) * n_smp)

do_test = function(n_ins, gene_len, rate)
    poisson.test(n_ins, gene_len * n_smp, rate)$p.value
result = tbl_df(as.data.frame(genes)) %>%
    mutate(p.value = purrr::map2_dbl(cis, width, do_test, rate=rate*2),
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

save(result, file=args$outfile)
