library(dplyr)
io = import('io')

# load sample-level aneuploidy w/ cutoffs
dset = io$load('../overview/aneuploidy_mad2.RData')
cis = io$read_table('171212 CIS RNA seq samples.tsv', header=TRUE) %>%
    mutate(Gene = stringr::str_trim(as.character(Gene))) %>%
    group_by(Sample, Gene) %>%
    summarize(n_ins = n()) %>%
    ungroup() %>%
    narray::construct(n_ins ~ Gene + Sample, fill=0) %>%
    narray::melt()
colnames(cis) = c("gene", "sample", "n_ins")

both = inner_join(dset, cis, by="sample") %>%
    select(-aneup, -mad2)

# FET of aneuploidy low/high vs. gene insertions
test_gene = function(df) {
    mat = data.matrix(df[,c('ins_gene','ins_total')])
    fisher.test(mat) %>%
        broom::tidy() %>%
        transmute(condition = df$condition[1],
                  log_odds = estimate,
                  p.value = p.value)
}

mad2 = filter(both, aneup_class != "high") %>%
    group_by(mad2_class, gene) %>%
    summarize(ins_gene = sum(n_ins)) %>%
    ungroup() %>%
    group_by(mad2_class) %>%
    mutate(ins_total = sum(ins_gene),
           condition = paste0("mad2-", mad2_class))
keep = mad2 %>%
    group_by(gene) %>%
    summarize(ins_gene = sum(ins_gene)) %>%
    filter(ins_gene >= 3) %>%
    pull(gene) %>%
    as.character()
mad2 = mad2[mad2$gene %in% keep,]

re = mad2 %>%
    group_by(gene) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, test_gene)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(p.value) %>%
    filter(adj.p < 0.2)
