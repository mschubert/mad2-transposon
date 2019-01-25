library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('d', 'dset', 'gene expression', '../expr_diff/eset_MILE.RData'),
    opt('s', 'geneset', 'gene set RData', '../data/genesets/human/GO_Biological_Process_2018.RData'),
    opt('o', 'outfile', 'save RData', 'gsva_mad2pb/GO_Biological_Process_2018.RData')
)

sets = io$load(args$geneset)
dset = io$load(args$dset)
scores = GSVA::gsva(dset$expr, sets, parallel.sz=1)
save(scores, file=args$outfile)
