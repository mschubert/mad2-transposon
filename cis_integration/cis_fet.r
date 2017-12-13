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
                  ins_cond = df$ins_gene[1],
                  ins_rest = df$ins_gene[2],
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
    filter(ins_gene >= 5) %>%
    pull(gene) %>%
    as.character()
mad2 = mad2[mad2$gene %in% keep,]

mad2 = mad2 %>%
    group_by(gene) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, test_gene)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(p.value) %>%
    filter(adj.p < 0.5)



aneup = filter(both, mad2_class == "low") %>%
    group_by(aneup_class, gene) %>%
    summarize(ins_gene = sum(n_ins)) %>%
    ungroup() %>%
    group_by(aneup_class) %>%
    mutate(ins_total = sum(ins_gene),
           condition = paste0("aneup-", aneup_class))
keep = aneup %>%
    group_by(gene) %>%
    summarize(ins_gene = sum(ins_gene)) %>%
    filter(ins_gene >= 5) %>%
    pull(gene) %>%
    as.character()
aneup = aneup[aneup$gene %in% keep,]

aneup = aneup %>%
    group_by(gene) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, test_gene)) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(p.value) %>%
    filter(adj.p < 0.5)

save(mad2, aneup, file="cis_fet.RData")
