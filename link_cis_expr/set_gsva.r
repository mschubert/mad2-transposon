library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression', '../data/rnaseq/assemble.RData'),
    opt('s', 'geneset', 'gene set RData', '../data/genesets/GO_Biological_Process_2018.RData'),
    opt('o', 'outfile', 'save RData', 'gsva_mad2pb/GO_Biological_Process_2018.RData')
)

expr = io$load(args$expr)
sets = io$load(args$geneset)

scores = GSVA::gsva(expr, sets, parallel.sz=1)

save(scores, file=args$outfile)
