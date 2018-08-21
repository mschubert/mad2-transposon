library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'all samples .RData', '../data/cis/cis_per_tumor.RData'),
    opt('r', 'reads', 'min reads to consider ins', '20'),
    opt('s', 'sheet', 'exclude/replace yaml', 'analysis_set.yaml'),
    opt('o', 'outfile', 'filtered samples & positions', 'analysis_set.RData'))

sheet = io$read_yaml(args$sheet)

# s&t different tumors: 184, 443 (RNA also both)
ins = io$load(args$infile) %>%
    dplyr::select(sample, chr, position, reads) %>% # ignore strand
    distinct() %>%
    filter(!sample %in% names(sheet$exclude)) %>%
    mutate(sample = ifelse(sample %in% names(sheet$replace),
                           unlist(sheet$replace)[sample], sample)) %>%
    filter(grepl("[0-9]{3}[st]", sample),
           reads >= as.integer(args$reads)) %>%
    arrange(sample, chr, position)

save(ins, file=args$outfile)
