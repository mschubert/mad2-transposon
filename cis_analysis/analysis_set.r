library(dplyr)
io = import('io')
sys = import('sys')

#' Local hopping correction by maximum number of reads in a window
#'
#' @param pos  Integer vector with positions
#' @param reads  number of reads
#' @param wsize  window size
#' @return  TRUE/FALSE if position is considered a local hop
is_local_hop = function(pos, reads, wsize) {
    cur = seq_along(pos)
    before = c(1, seq_along(pos)[-length(pos)])
    after = c(seq_along(pos)[-1], length(pos))
    (reads[cur] < reads[before] & abs(pos[cur] - pos[before]) < wsize) |
        (reads[cur] < reads[after] & abs(pos[cur] - pos[after]) < wsize)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'all samples .RData', '../data/cis/cis_per_tumor.RData'),
    opt('r', 'reads', 'min reads to consider ins', '20'),
    opt('f', 'read_frac', 'min frac of reads for gene', '0.01'),
    opt('s', 'sheet', 'exclude/replace yaml', 'analysis_set.yaml'),
    opt('w', 'wsize', 'window size local hops', '0'),
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
    group_by(sample) %>%
    filter(reads >= sum(reads) * as.numeric(args$read_frac)) %>%
    ungroup() %>%
    arrange(sample, chr, position) %>%
    group_by(sample, chr) %>%
    mutate(is_local = is_local_hop(position, reads, as.integer(args$wsize))) %>%
    ungroup()

save(ins, file=args$outfile)
