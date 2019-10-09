io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'other analysis input', '../data/genesets/mouse/CH.HALLMARK.RData'),
    opt('o', 'outfile', 'rds', 'sig_myc.rds'))

dset = io$load(args$infile)$HALLMARK_MYC_TARGETS_V2

saveRDS(dset, file=args$outfile)
