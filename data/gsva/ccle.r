io = import('io')
sys = import('sys')
ccle = import('data/ccle')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/CH.HALLMARK.RData'),
    opt('t', 'threads', 'num', '12'),
    opt('o', 'outfile', 'rds', 'ccle/CH.HALLMARK.rds')
)

sets = io$load(args$setfile)
expr = ccle$basal_expression()

scores = GSVA::gsva(expr, sets, parallel.sz=as.integer(args$threads))

saveRDS(scores, file=args$outfile)
