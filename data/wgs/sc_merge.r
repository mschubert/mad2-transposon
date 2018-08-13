io = import('io')
sys = import('sys')
aneuf = import("tools/aneufinder")

args = sys$cmd$parse(
    opt('o', 'outfile', 'merged .RData', 'sc_merge.RData'),
    arg('infiles', 'scWGS .RData', sprintf("%s.RData", c("T401", "T419", "S413")), arity='*'))

merged = io$load(args$infiles) %>%
    lapply(aneuf$consensus_ploidy) %>%
    dplyr::bind_rows(.id = "sample")

save(merged, file=args$outfile)
