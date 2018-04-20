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
    opt('s', 'samples', 'min samples w/ insertions in gene', '5'))

cis = io$load("dset.RData")
samples = cis %>%
    select(sample, aneup) %>%
    distinct() %>%
    mutate(reads = 0)

do_fit = function(data) {
    no_reads = samples %>% filter(!sample %in% data$sample)
    data2 = dplyr::bind_rows(data, no_reads)
    stopifnot(nrow(data2) == nrow(samples))

    betareg(aneup ~ reads, data=data2) %>%
        broom::tidy() %>%
        filter(term == "reads") %>%
        select(-component, -term)
}

assoc = . %>%
    group_by(sample) %>%
    mutate(reads = reads / sum(reads)) %>%
    group_by(sample, gene_name, ensembl_gene_id, aneup) %>%
    summarize(reads = sum(reads)) %>%
    group_by(gene_name, ensembl_gene_id) %>%
    filter(sum(!is.na(aneup)) >= as.integer(args$sample)) %>%
    tidyr::nest() %>%
    mutate(fit = purrr::map(data, do_fit)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    arrange(p.value) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

hits = assoc(filter(cis, type == "hit"))
hits_cancer = assoc(filter(cis, type == "hit", known_cancer))
near = assoc(cis)
near_cancer = assoc(filter(cis, known_cancer))

save(hits, hits_cancer, near, near_cancer, file="betareg.RData")
# plot volcano, top hit fits
