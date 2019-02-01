library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('g', 'genes', '.RData', 'genes_mad2pb.RData'),
    opt('s', 'sets', '.RData', 'sets_mad2pb.RData'),
    opt('v', 'tfs', '.RData', 'tfs_mad2pb.RData'),
    opt('m', 'merge', '.yaml', 'stat_merge.yaml'),
    opt('o', 'outfile', 'RData', 'both_mad2pb.RData')
)

mfile = io$read_yaml(args$merge)

genes = io$load(args$genes)
sets = io$load(args$sets)
tfs = io$load(args$tfs)

stopifnot(Reduce(all.equal, list(genes$groups, sets$groups, tfs$groups)))
expr = rbind(genes$expr[mfile$genes,, drop=FALSE],
             sets$expr[mfile$sets, colnames(genes$expr), drop=FALSE],
             tfs$expr[mfile$tfs, colnames(genes$expr), drop=FALSE])
groups = genes$groups

save(expr, groups, file=args$outfile)
