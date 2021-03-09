library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'dset.rds'),
    opt('g', 'go', 'rds', '../data/gsva/tcga-brca/GO_Biological_Process_2020.rds'),
    opt('p', 'plotfile', 'pdf', 'stat1ko_cor.pdf')
)

ds = readRDS(args$dset)
dset = cbind(ds$meta[c("purity")], as.data.frame(ds$dmat))
go = t(readRDS(args$go))
narray::intersect(dset, go, along=1)

test_one = function(cn, mat) {
    lm(mat[,cn] ~ rev24_stat1_over_wt, data=dset) %>%
        broom::tidy() %>%
        filter(term != "(Intercept)")
}
res_go = sapply(colnames(go), test_one, mat=go, simplify=FALSE) %>%
    bind_rows(.id="label") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

res_other = setdiff(colnames(dset), "rev24_stat1_over_wt") %>%
    sapply(test_one, mat=dset, simplify=FALSE) %>%
    bind_rows(.id="label") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pdf(args$plotfile, 10, 8)
res_go %>%
    mutate(size = 1) %>%
    plt$p_effect("adj.p", thresh=1e-10) %>%
    plt$volcano(text.size=2.5, repel=TRUE, label_top=50, p=1e-10, pos_label_bias=0.5, max.overlaps=50) +
        ggtitle("GO bio bp 2020 <> rev24_stat1_over_wt")

res_other %>%
    mutate(size = 1) %>%
    plt$p_effect("adj.p", thresh=1e-10) %>% # 1e-4: macrophages huge effect size, but plot doesn;t look llike
    plt$volcano(text.size=2.5, repel=TRUE, label_top=50, p=1e-10, pos_label_bias=0.5, max.overlaps=50) +
        ggtitle("rev24_stat1_over_wt")
dev.off()
