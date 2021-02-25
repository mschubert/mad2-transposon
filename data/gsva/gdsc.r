io = import('io')
sys = import('sys')
gdsc = import('data/gdsc')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/CH.HALLMARK.RData'),
    opt('t', 'threads', 'num', '12'),
    opt('o', 'outfile', 'rds', 'gdsc/CH.HALLMARK.rds')
)

sets = io$load(args$setfile)
expr = gdsc$basal_expression()

scores = GSVA::gsva(expr, sets, parallel.sz=as.integer(args$threads))

saveRDS(scores, file=args$outfile)
