io = import('io')
sys = import('sys')
gdsc = import('data/gdsc')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/CH.HALLMARK.RData'),
    opt('o', 'outfile', 'rds', 'gdsc/CH.HALLMARK.rds'))

sets = io$load(args$setfile)
expr = gdsc$basal_expression()

scores = GSVA::gsva(expr, sets, parallel.sz=1)

saveRDS(scores, file=args$outfile)
