library(dplyr)
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dna', 'copy number segments', '../ploidy_from_wgs/copy_segments.RData'),
    opt('f', 'fractions', 'fraction merge tsv', 'wgs_merge.tsv'),
#    opt('r', 'rna', 'rna copy segments', '../ploidy_from_rnaseq/eT_ploidy.RData'),
#    opt('s', 'scdna', 'single cell measures', '../data/wgs/'),
    opt('o', 'outfile', 'results RData', 'gene_copies.RData'))

fracs = io$read_table(args$fractions, header=TRUE) %>%
    select(sample, subset, weight)

coords = seq$coords$gene("ensembl_gene_id", granges=TRUE,
    dset="mmusculus_gene_ensembl", chromosomes=c(1:19,'X'))

gene_copies = readRDS(args$dna)$segments %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    seq$intersect(coords) %>%
    select(subset=sample, ensembl_gene_id, ploidy) %>%
    mutate(sample = sub("(-high)|(-low)", "", subset)) %>%
    left_join(fracs) %>%
    mutate(weight = ifelse(sample == subset, 1, weight)) %>%
    group_by(sample, ensembl_gene_id) %>%
    summarize(ploidy = weighted.mean(ploidy, weight)) %>%
    ungroup() %>%
    na.omit() %>% # drop 417s, 446s[nofrac], 449s, 462s, 477t, 489s, 613t
    narray::construct(ploidy ~ ensembl_gene_id + sample)

#rna = readRDS(args$rna)$segments %>%
#    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
#    seq$intersect(coords) %>%
#    select(sample = sample, ensembl_gene_id, ploidy=expr)

saveRDS(gene_copies, file=args$outfile)
