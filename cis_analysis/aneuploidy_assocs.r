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

#' KS statistic for aneuploidy associations
ks_test = function() {
}

args = sys$cmd$parse(
    opt('a', 'aneup', 'aneuploidy .tsv', '../ploidy_compare/compare_ploidy.tsv'),
    opt('p', 'poisson', 'cis assocs RData', 'poisson.RData'),
    opt('o', 'outfile', 'aneuploidy assocs', 'aneuploidy_assocs.RData'))

dset = io$load(args$poisson)
genes = dset$result %>%
    filter(adj.p < 1e-5) %>%
    pull(external_gene_name)

# load and merge aneuploidy scores
aneup = io$read_table(args$aneup, header=TRUE) %>%
    select(-coverage, -tissue) %>%
    filter(!duplicated(data.frame(type, sample)),
           !sample %in% c("S", "T")) %>%
    mutate(sample = fix_ids(sample)) %>%
    filter(sample %in% aset$sample) %>%
    tidyr::spread("type", "aneuploidy") %>%
    group_by(sample) %>%
    summarize(aneup = `WGS (merged)` %or% `WGS (30-cell)` %or% `RNA-seq (eT)`) %>%
    arrange(-aneup)


