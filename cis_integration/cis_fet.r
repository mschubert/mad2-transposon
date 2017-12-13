library(dplyr)
library(cowplot)
io = import('io')

# load sample-level aneuploidy w/ cutoffs
dset = io$load('aneuploidy_mad2.RData')
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
           condition = paste0("mad2-", mad2_class)) %>%
    ungroup()
mad2_ins = mad2 %>%
    select(condition, ins_total) %>%
    distinct()
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
    arrange(p.value)

p1 = ggplot(mad2_ins, aes(x=condition, y=ins_total, fill=condition)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    guides(fill=FALSE)
p2 = mad2 %>%
    head(10) %>%
    mutate(gene = factor(gene, levels=unique(gene))) %>%
    tidyr::gather("cond", "n_ins", ins_cond, ins_rest)
p2$condition[p2$cond == "ins_rest"] = sub("high", "low", p2$condition[p2$cond=="ins_cond"])
p2 = ggplot(p2, aes(x=gene, y=n_ins, fill=condition)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
p_mad2 = plot_grid(p1, p2, ncol=2, rel_widths=c(1,3))


aneup = filter(both, mad2_class == "low") %>%
    group_by(aneup_class, gene) %>%
    summarize(ins_gene = sum(n_ins)) %>%
    ungroup() %>%
    group_by(aneup_class) %>%
    mutate(ins_total = sum(ins_gene),
           condition = paste0("aneup-", aneup_class)) %>%
    ungroup()
aneup_ins = aneup %>%
    select(condition, ins_total) %>%
    distinct()
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
    arrange(p.value)

p1 = ggplot(aneup_ins, aes(x=condition, y=ins_total, fill=condition)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    guides(fill=FALSE)
p2 = aneup %>%
    head(10) %>%
    mutate(gene = factor(gene, levels=unique(gene))) %>%
    tidyr::gather("cond", "n_ins", ins_cond, ins_rest)
p2$condition[p2$cond == "ins_rest"] = sub("high", "low", p2$condition[p2$cond=="ins_cond"])
p2$n_ins[p2$gene == "Intergenic"] = p2$n_ins[p2$gene == "Intergenic"] / 50
levels(p2$gene)[levels(p2$gene) == "Intergenic"] = "Intergenic / 50"
p2 = ggplot(p2, aes(x=gene, y=n_ins, fill=condition)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
p_aneup = plot_grid(p1, p2, ncol=2, rel_widths=c(1,3))

pdf("cis_fet.pdf")
print(p_mad2)
print(p_aneup)
dev.off()


save(mad2, mad2_ins, aneup, aneup_ins, file="cis_fet.RData")
