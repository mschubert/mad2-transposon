library(dplyr)
library(tidygraph)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
io = import('io')

genes = c("Ets1", "Erg", "Stat1", "Pias1")
dset = io$load("../expr_diff/eset_Mad2PB.RData")

meta = as.data.frame(SummarizedExperiment::colData(dset$eset))

expr = reshape2::melt(dset$vs[genes,]) %>%
    dplyr::rename(gene=Var1, sample=Var2, expr=value) %>%
    inner_join(meta %>% select(sample, type)) %>%
    filter(type != "unknown") %>%
    mutate(type = as.character(type))

eplot = ggplot(expr, aes(x=type, y=expr)) +
    geom_boxplot(aes(fill=type), outlier.shape=NA) +
    facet_wrap(~ gene, scales="free_y", ncol=1)

pdf("Fig3-RNA+EtsErg.pdf", 4, 3.5)
#FIXME: annoying empty plots
c1 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "Myeloid"])), main="Myeloid")
c2 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "Other"])), main="B-like")
c3 = corrplot::corrplot(cor(t(dset$vs[genes, meta$type == "T-cell"])), main="T-ALL")
#{ c1 + c2 + c3 + plot_layout(ncol=1) } + { eplot } + plot_layout(nrow=1)
dev.off()

pdf("Fig3-RNA+EtsErg_2.pdf", 4, 10)
print(eplot + labs(x="Type", y="Gene expression"))
dev.off()
