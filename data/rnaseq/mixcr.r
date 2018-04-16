io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'tsv to save to', 'mixcr_Mad2+MB.tsv'),
    arg('infiles', 'mixcr result files', arity='*'))

combined = args$infiles %>%
    lapply(io$read_table, header=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(args$infiles))) %>%
    dplyr::bind_rows(.id="sample") %>%
    dplyr::mutate(type = substr(allVHitsWithScore, 1, 4))

write.table(combined, sep="\t", file=args$outfile, quote=FALSE, row.names=FALSE)
