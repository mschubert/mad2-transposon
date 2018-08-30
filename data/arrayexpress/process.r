library(dplyr)
ma = import('process/microarray')
idmap = import('process/idmap')
sys = import('sys')

args = sys$cmd$parse(
    opt('a', 'accession', 'E-GEOD/MTAB/etc.', 'E-GEOD-28497'),
    opt('s', 'summarize', 'identifier type', 'ensembl_gene_id'),
    opt('o', 'outfile', 'save to .RData', 'E-GEOD-28497.RData'))

expr = ArrayExpress::ArrayExpress(args$accession) %>%
    ma$qc() %>%
    ma$normalize() %>%
    ma$annotate(summarize=args$summarize)

save(expr, file=args$outfile)
