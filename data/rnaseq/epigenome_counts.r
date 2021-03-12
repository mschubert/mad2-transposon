library(dplyr)
library(plyranges)
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'epigenome', 'tsv', 'epigenome_ensembl102.tsv'),
    opt('o', 'outfile', 'rds', 'epigenome_counts.rds'),
    arg('bamfiles', 'bam', arity='*',
        list.files("aligned", "\\.bam$", recursive=TRUE, full.names=TRUE))
)

epi = readr::read_tsv(args$epigenome) %>%
    select(-`Epigenome name`, -Activity) %>%
    distinct() %>%
    makeGRangesFromDataFrame(seqnames.field="Chromosome/scaffold name",
                             start.field="Start (bp)", end.field="End (bp)",
                             keep.extra.columns=TRUE)

reads = matrix(NA, ncol=length(args$bamfiles), nrow=length(epi),
               dimnames=list(epi$`Regulatory stable ID`, args$bamfiles))

for (bam in args$bamfiles) {
    gr = seq$read_granges(bam)
    reads[,bam] = count_overlaps(epi, gr)
}

saveRDS(list(epi=epi, reads=reads), file=args$outfile)
