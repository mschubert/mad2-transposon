# logistic regression
#
# 20(ish) read cutoff per transposon insertion
# aneuploidy scores as continuous variable
# test whether or not each gene has insertion w/ increasing aneuploidy score
library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('r', 'reads', 'min number of reads', '20'),
    opt('s', 'samples', 'min samples w/ insertions in gene', '5'))

cis = io$load("dset.RData")
samples = cis %>%
    select(sample, aneup, n_ins_smp) %>%
    distinct() %>%
    mutate(reads = 0)

do_fit = function(data) {
    no_reads = samples %>% filter(!sample %in% data$sample)
    data2 = dplyr::bind_rows(data, no_reads) %>%
        mutate(insertion = reads >= as.integer(args$reads))
    stopifnot(nrow(data2) == nrow(samples))

    glm(insertion ~ n_ins_smp + aneup, family=binomial(link='logit'), data=data2) %>%
        broom::tidy() %>%
        filter(term == "aneup") %>%
        select(-term)
}

assoc = . %>%
    group_by(sample, gene_name, ensembl_gene_id, aneup, n_ins_smp) %>%
    summarize(reads = sum(reads)) %>%
    group_by(gene_name, ensembl_gene_id) %>%
    filter(sum(!is.na(aneup) & reads >= as.integer(args$reads)) >=
           as.integer(args$sample)) %>%
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

save(hits, hits_cancer, near, near_cancer, file="logit.RData")
