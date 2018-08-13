library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')

#' source sample IDs are both 123S and T567
fix_ids = function(x) {
    letter = tolower(sub("[0-9]+", "", x))
    nums = sub("[ST]", "", x)
    paste0(nums, letter)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'aneuploidy .tsv', 'compare_ploidy.tsv'),
    opt('m', 'meta', 'sample types', '../data/meta/meta.RData'),
    opt('o', 'outfile', 'aneuploidy assocs', 'analysis_set.RData'))

meta = io$load(args$meta)$meta %>%
    transmute(sample = `Hist nr.`,
              tissue = rep(list(c("spleen", "thymus")), nrow(.)),
              type = ifelse(is.na(Diagnosis), "NA", Diagnosis),
              type = sub("(^[^ ]+).*", "\\1", type)) %>%
    tidyr::unnest() %>%
    mutate(sample = paste0(sample, substr(tissue, 1, 1)))

#TODO: track which data source we use to determine aneuploidy
aneup = io$read_table(args$infile, header=TRUE) %>%
    select(-coverage, -tissue) %>%
    filter(!duplicated(data.frame(type, sample)),
           !sample %in% c("S", "T")) %>%
    mutate(sample = fix_ids(sample)) %>%
    tidyr::spread("type", "aneuploidy") %>%
    group_by(sample) %>%
    summarize(aneup = `WGS (merged)` %or% `WGS (30-cell)` %or% `RNA-seq (eT)`) %>%
    arrange(-aneup) %>%
    left_join(meta)

save(aneup, file=args$outfile)
