sys = import('sys')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/MSigDB_Hallmark_2020.rds'),
    opt('t', 'threads', 'num', '12'),
    opt('o', 'outfile', 'rds', 'mad2pb/MSigDB_Hallmark_2020.rds')
)

sets = readRDS(args$setfile)
expr = readRDS("../../expr_diff/eset_Mad2PB.rds")$vs

scores = GSVA::gsva(expr, sets, parallel.sz=as.integer(args$threads))

saveRDS(scores, file=args$outfile)
