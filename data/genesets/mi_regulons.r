io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('n', 'network', 'aracne results RData', '../networks/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'gene set RData', 'MI_regulons.RData'))

regs = unstack(io$load(args$network)[c(2,1)])

save(regs, file=args$outfile)
