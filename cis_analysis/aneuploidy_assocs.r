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
ks_test = function(data) {
    mod = ks.test(data$aneup[data$reads], data$aneup[!data$reads])
    mod$model = select(data, sample, aneup, reads)
    mod %>%
        broom::tidy() %>%
        mutate(estimate = mean(data$aneup[data$reads]) - mean(data$aneup),
               statistic = sign(estimate) / p.value,
               size = sum(data$reads),
               mod = list(mod))
}

args = sys$cmd$parse(
    opt('a', 'aneup', 'aneuploidy .tsv', '../ploidy_compare/analysis_set.RData'),
    opt('p', 'poisson', 'cis assocs RData', 'poisson.RData'),
    opt('o', 'outfile', 'aneuploidy assocs', 'aneuploidy_assocs.RData'))

aneup = io$load(args$aneup)
dset = io$load(args$poisson)
genes = dset$result %>%
    filter(adj.p < 1e-5) %>%
    pull(external_gene_name)

result = dset$samples %>%
    select(-n_ins) %>%
    filter(external_gene_name %in% genes) %>%
    mutate(reads=TRUE) %>%
    tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
    inner_join(aneup, by="sample") %>%
    group_by(external_gene_name) %>%
    tidyr::nest() %>%
    mutate(fit = purrr::map(data, ks_test)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    select(external_gene_name, size, estimate, statistic, p.value, adj.p) %>%
    arrange(adj.p, p.value)

save(result, file=args$outfile)