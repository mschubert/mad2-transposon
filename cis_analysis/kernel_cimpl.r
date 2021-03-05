library(dplyr)
library(cimpl)
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'cis rds', '../data/cis/cis_per_tumor.rds'),
    opt('u', 'upstream', 'bp to include before gene', '10000'),
    opt('d', 'downstream', 'bp to include after gene', '0'),
    opt('r', 'reads', 'min num of reads', '20'),
    opt('o', 'outfile', 'insertion statistics rds', 'kernel_cimpl.rds')
)

ins = io$load(args$infile) %>%
    dplyr::select(sample, chr, position, reads) %>% # ignore strand
    distinct() %>%
    filter(reads >= as.integer(args$reads)) %>%
    transmute(sampleID = sample,
              location = position,
              contig_depth = reads,
              chr = paste0("chr", chr))

run = doCimplAnalysis(ins, BSgenome=BSgenome.Mmusculus.UCSC.mm10::Mmusculus,
                      system="PB", lhc.method="exclude", verbose=TRUE)
cis = getCISs(run, order.by="p_value", decreasing=FALSE)

#TODO: map to genes, subset similarly to poisson
