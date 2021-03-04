library(dplyr)
sys = import('sys')
aneuf = import("tools/aneufinder")

args = sys$cmd$parse(
    opt('o', 'outfile', 'merged .RData', 'sc_merge.rds'),
    arg('infiles', 'scWGS rds', arity='*',
        sub(".yaml", ".rds", list.files(pattern="[0-9]{3}[st].yaml$"), fixed=TRUE))
)

merged = lapply(args$infiles, readRDS) %>%
    lapply(aneuf$consensus_ploidy) %>%
    dplyr::bind_rows(.id = "sample")

saveRDS(merged, file=args$outfile)
