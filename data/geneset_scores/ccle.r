io = import('io')
sys = import('sys')
ccle = import('data/ccle')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/CH.HALLMARK.RData'),
    opt('o', 'outfile', 'rds', 'ccle/CH.HALLMARK.rds'))

sets = io$load(args$setfile)
expr = ccle$basal_expression()

scores = GSVA::gsva(expr, sets, parallel.sz=1)

saveRDS(scores, file=args$outfile)
