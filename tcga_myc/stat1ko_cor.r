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

test_one = function(cn) {
    lm(go[,cn] ~ rev24_stat1_over_wt, data=dset) %>%
        broom::tidy() %>%
        filter(term != "(Intercept)")
}
res = sapply(colnames(go), test_one, simplify=FALSE) %>%
    bind_rows(.id="label") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pdf(args$plotfile, 10, 8)
res %>%
    mutate(size = 1) %>%
    plt$p_effect("adj.p", thresh=0.01) %>%
    plt$volcano(text.size=2.5, repel=TRUE, label_top=50, pos_label_bias=0.5, max.overlaps=50)
dev.off()
