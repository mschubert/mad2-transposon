library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'gene expression', '../expr_diff/eset_Mad2PB.RData'),
    opt('s', 'geneset', 'gene set RData', '../data/genesets/mouse/GO_Biological_Process_2018.RData'),
    opt('o', 'outfile', 'save RData', 'gsva_mad2pb/GO_Biological_Process_2018.RData')
)

dset = io$load(args$dset)
sets = io$load(args$geneset)
scores = GSVA::gsva(dset$vs, sets, parallel.sz=1)
save(scores, file=args$outfile)
