library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('g', 'genes', '.RData', 'genes_mad2pb.RData'),
    opt('s', 'sets', '.RData', 'sets_mad2pb.RData'),
    opt('o', 'outfile', 'RData', 'both_mad2pb.RData')
)

genes = io$load(args$genes)
sets = io$load(args$sets)

stopifnot(all.equal(genes$groups, sets$groups))
expr = rbind(genes$expr, sets$expr[,colnames(genes$expr)])
groups = genes$groups

save(expr, groups, file=args$outfile)
