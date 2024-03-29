# logistic regression
#
# 20(ish) read cutoff per transposon insertion
# aneuploidy scores as continuous variable
# test whether or not each gene has insertion w/ increasing aneuploidy score
library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('r', 'reads', 'min number of reads', '20'),
    opt('s', 'samples', 'min samples w/ insertions in gene', '5'),
    opt('o', 'outfile', 'file to save to', 'gamma_adjNins.RData'))

cis = io$load("dset.RData")
samples = cis %>%
    select(sample, aneup, n_ins_smp) %>%
    distinct() %>%
    mutate(reads = 0)

do_fit = function(data) {
    no_reads = samples %>% filter(!sample %in% data$sample)
    data2 = dplyr::bind_rows(data, no_reads) %>%
        mutate(reads = reads >= as.integer(args$reads))
    stopifnot(nrow(data2) == nrow(samples))

    mod = glm(aneup ~ n_ins_smp + reads, family=Gamma(), data=data2) %catch%
        return(data.frame(estimate=NA))
    mod$model = na.omit(select(data2, sample, aneup, reads))
    mod %>%
        broom::tidy() %>%
        filter(term == "readsTRUE") %>%
        select(-term) %>%
        mutate(size = length(unique(data$sample)),
               mod = list(mod))
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
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    select(-mod, mod)

hits = assoc(filter(cis, hit_dist == 0))
hits_cancer = assoc(filter(cis, hit_dist == 0, known_cancer))
near = assoc(cis)
near_cancer = assoc(filter(cis, known_cancer))

save(hits, hits_cancer, near, near_cancer, file=args$outfile)
