sys = import('sys')
ccle = import('data/ccle')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/MSigDB_Hallmark_2020.rds'),
    opt('t', 'threads', 'num', '12'),
    opt('o', 'outfile', 'rds', 'ccle/MSigDB_Hallmark_2020.rds')
)

sets = readRDS(args$setfile)
expr = ccle$basal_expression()

scores = GSVA::gsva(expr, sets, parallel.sz=as.integer(args$threads))

saveRDS(scores, file=args$outfile)
