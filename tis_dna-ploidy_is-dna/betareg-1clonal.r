# beta regression
#
# use the fraction of reads that are in transposon
# aneuploidy scores as continuous variable
# test whether or not each gene has insertion w/ increasing aneuploidy score
library(dplyr)
library(betareg)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('s', 'samples', 'min samples w/ insertions in gene', '5'),
    opt('d', 'decay', 'decay factor for flanking genes', '0.7'),
    opt('o', 'outfile', 'file to save to', 'betareg-1clonal.RData'))

cis = io$load("dset.RData") %>%
    group_by(sample) %>%
    mutate(reads = as.numeric(args$decay)^hit_dist * reads / max(reads))

samples = cis %>%
    select(sample, aneup) %>%
    distinct() %>%
    mutate(reads = 0)

do_fit = function(data) {
    no_reads = samples %>% filter(!sample %in% data$sample)
    data2 = dplyr::bind_rows(data, no_reads)
    stopifnot(nrow(data2) == nrow(samples))

    mod = betareg(aneup ~ reads, data=data2)
    mod %>%
        broom::tidy() %>%
        filter(term == "reads") %>%
        select(-component, -term) %>%
        mutate(size = length(unique(data$sample)),
               mod = list(mod))
}

assoc = . %>%
    group_by(sample, gene_name, ensembl_gene_id, aneup) %>%
    summarize(reads = sum(reads)) %>%
    group_by(gene_name, ensembl_gene_id) %>%
    filter(sum(!is.na(aneup)) >= as.integer(args$sample)) %>%
    tidyr::nest() %>%
    mutate(fit = purrr::map(data, do_fit)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    arrange(p.value) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    select(-mod, mod)

hits = assoc(filter(cis, hit_dist == 0))
hits_cancer = assoc(filter(cis, hit_dist == 0, known_cancer))
near = assoc(cis)
near_cancer = assoc(filter(cis, known_cancer))

save(hits, hits_cancer, near, near_cancer, file=args$outfile)
